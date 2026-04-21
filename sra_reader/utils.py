"""
Utilities: dependency checks, SRA metadata extraction, vdb-dump launcher.
"""

import csv
import io
import subprocess
import shutil
import re
from dataclasses import dataclass
from typing import Optional


class SRAReaderError(Exception):
    """Base exception for sra_reader."""
    pass


class DependencyError(SRAReaderError):
    """Raised when required CLI tools are missing."""
    pass


class SRAFormatError(SRAReaderError):
    """Raised when SRA file is malformed or unreadable."""
    pass


REQUIRED_TOOLS = ["vdb-dump", "sra-stat"]
OPTIONAL_TOOLS = ["prefetch"]


def check_dependencies(raise_on_missing: bool = True) -> dict[str, Optional[str]]:
    """
    Check that required SRA toolkit binaries are on PATH.

    Returns dict of tool_name -> path (or None if missing).
    Raises DependencyError if raise_on_missing and any required tool is absent.
    """
    found = {}
    missing = []
    for tool in REQUIRED_TOOLS + OPTIONAL_TOOLS:
        path = shutil.which(tool)
        found[tool] = path
        if path is None and tool in REQUIRED_TOOLS:
            missing.append(tool)
    if missing and raise_on_missing:
        raise DependencyError(
            f"Missing required SRA toolkit binaries: {', '.join(missing)}. "
            f"Install via: conda install -c bioconda sra-tools"
        )
    return found


@dataclass
class SRAInfo:
    """Metadata about an SRA accession or file."""
    accession: str
    num_spots: int
    num_reads_per_spot: int
    read_lengths: list[int]
    read_types: list[str]
    platform: str
    is_variable_length: bool

    @property
    def num_segments(self) -> int:
        return self.num_reads_per_spot


def sra_info(sra_path: str) -> SRAInfo:
    """
    Extract structural metadata from an SRA file/accession.

    Uses sra-stat for spot count and vdb-dump for read structure.
    """
    check_dependencies()
    num_spots = _get_spot_count(sra_path)
    platform = _get_platform(sra_path)
    read_types, read_lengths, is_variable = _get_read_structure(sra_path)
    accession = sra_path.rstrip("/").split("/")[-1].replace(".sra", "")

    return SRAInfo(
        accession=accession,
        num_spots=num_spots,
        num_reads_per_spot=len(read_types),
        read_lengths=read_lengths,
        read_types=read_types,
        platform=platform,
        is_variable_length=is_variable,
    )


def _get_spot_count(sra_path: str) -> int:
    """Get total spot count using sra-stat."""
    try:
        result = subprocess.run(
            ["sra-stat", "--xml", "--quick", sra_path],
            capture_output=True, text=True, timeout=60,
        )
        if result.returncode != 0:
            raise SRAFormatError(
                f"sra-stat failed on {sra_path}: {result.stderr.strip()}"
            )
        match = re.search(r'spot_count="(\d+)"', result.stdout)
        if match:
            return int(match.group(1))
        raise SRAFormatError("Could not parse spot_count from sra-stat output")
    except subprocess.TimeoutExpired:
        raise SRAReaderError(f"sra-stat timed out on {sra_path}")
    except FileNotFoundError:
        raise DependencyError("sra-stat not found")


def _get_platform(sra_path: str) -> str:
    """Extract platform from vdb-dump --info."""
    try:
        result = subprocess.run(
            ["vdb-dump", "--info", sra_path],
            capture_output=True, text=True, timeout=30,
        )
        for line in result.stdout.splitlines():
            if "platf" in line.lower():
                return line.split(":")[-1].strip()
        return "UNKNOWN"
    except Exception:
        return "UNKNOWN"


def _parse_vdb_csv_line(line: str) -> list[str]:
    """
    Parse a vdb-dump CSV line respecting quoted fields.

    vdb-dump -f csv produces lines like:
      ACGTACGT,"SRA_READ_TYPE_BIOLOGICAL, SRA_READ_TYPE_BIOLOGICAL","150, 150","0, 150"
    """
    reader = csv.reader(io.StringIO(line))
    for row in reader:
        return row
    return []


def _parse_read_type(type_str: str) -> str:
    """
    Convert SRA read type string or int to 'Biological' or 'Technical'.
    """
    s = type_str.strip()
    if "BIOLOGICAL" in s.upper():
        return "Biological"
    if "TECHNICAL" in s.upper():
        return "Technical"
    try:
        val = int(s)
        return "Biological" if (val & 1) else "Technical"
    except ValueError:
        return s


def _get_read_structure(
    sra_path: str,
    sample_spots: int = 10,
) -> tuple[list[str], list[int], bool]:
    """
    Determine read structure by sampling first N spots.

    Returns (read_types, read_lengths, is_variable_length).
    """
    result = subprocess.run(
        [
            "vdb-dump", sra_path,
            "-C", "READ_TYPE,READ_LEN",
            "-R", f"1-{sample_spots}",
            "-f", "csv",
        ],
        capture_output=True, text=True, timeout=30,
    )
    if result.returncode != 0:
        raise SRAFormatError(
            f"vdb-dump failed reading structure: {result.stderr.strip()}"
        )

    all_types = []
    all_lengths = []

    for line in result.stdout.strip().splitlines():
        if not line:
            continue
        fields = _parse_vdb_csv_line(line)
        if len(fields) < 2:
            continue
        type_labels = [_parse_read_type(t) for t in fields[0].split(",")]
        lengths = [int(l.strip()) for l in fields[1].split(",")]
        all_types.append(type_labels)
        all_lengths.append(lengths)

    if not all_types:
        raise SRAFormatError(f"No spots found in {sra_path}")

    read_types = all_types[0]
    read_lengths = all_lengths[0]

    is_variable = False
    for lens in all_lengths[1:]:
        if lens != read_lengths:
            is_variable = True
            break

    return read_types, read_lengths, is_variable


def run_vdb_dump(
    sra_path: str,
    columns: list[str],
    row_range: Optional[tuple[int, int]] = None,
    fmt: str = "csv",
) -> subprocess.Popen:
    """
    Launch vdb-dump as a streaming subprocess.

    Returns a Popen object whose stdout can be iterated line-by-line.
    Uses binary mode (stdout=PIPE without text=True) for Cython fast path.
    """
    cmd = ["vdb-dump", sra_path, "-C", ",".join(columns), "-f", fmt]
    if row_range is not None:
        cmd.extend(["-R", f"{row_range[0]}-{row_range[1]}"])
    return subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=4 * 1024 * 1024,  # 4MB buffer
    )


def has_xsra() -> bool:
    """Check if xsra (Arc Institute's fast SRA reader) is available."""
    return shutil.which("xsra") is not None


def run_xsra_stream(
    sra_path: str,
    segment_indices: Optional[list[int]] = None,
    row_limit: Optional[int] = None,
    n_threads: int = 4,
    fmt: str = "fasta",
) -> subprocess.Popen:
    """
    Launch xsra as a streaming subprocess for reading SRA directly.

    xsra uses the VDB C library directly (no text serialization) and
    supports multi-threaded decompression. This is 5-20x faster than
    vdb-dump for streaming reads.

    Install: cargo install xsra
    (requires Rust toolchain)

    Output format:
      --fasta: >spot_id.segment_idx\\nSEQUENCE
      --fastq: @spot_id.segment_idx\\nSEQUENCE\\n+\\nQUALITY

    Args:
        sra_path: Path to .sra file.
        segment_indices: Which segments to extract (0-based). None = all.
        row_limit: Max spots to read. None = all.
        n_threads: Thread count for xsra (0 = auto).
        fmt: "fasta" or "fastq".

    Returns:
        Popen with stdout streaming the output.
    """
    cmd = ["xsra", "dump", sra_path, f"-T{n_threads}"]

    if fmt == "fasta":
        cmd.append("-fa")
    else:
        cmd.append("-fq")

    if segment_indices is not None:
        cmd.extend(["-I", ",".join(str(i) for i in segment_indices)])

    if row_limit is not None:
        cmd.extend(["-l", str(row_limit)])

    return subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=4 * 1024 * 1024,
    )

