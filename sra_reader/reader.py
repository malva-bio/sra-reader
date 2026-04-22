"""
Core SRA reader: streams reads directly from .sra files.

The hot loop uses the Cython _parser module for line parsing.
Structure is detected once at init; the fast path only fetches
the READ column and slices with precomputed offsets in C.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from typing import Iterator, Optional

from .utils import (
    SRAInfo,
    SRAReaderError,
    SRAFormatError,
    check_dependencies,
    run_vdb_dump,
    run_xsra_stream,
    has_xsra,
    sra_info,
    _parse_vdb_csv_line,
    _parse_read_type,
)
from ._parser import parse_line_fixed_seqonly, parse_line_fixed


@dataclass(frozen=True, slots=True)
class SRASegment:
    """A single read segment (R1/R2/I1/I2 etc.)."""
    sequence: str
    quality: str
    read_type: str
    segment_index: int


@dataclass(slots=True)
class SRASpot:
    """A spot containing all read segments."""
    spot_id: int
    segments: list[SRASegment] = field(default_factory=list)

    @property
    def num_segments(self) -> int:
        return len(self.segments)

    def biological_segments(self) -> list[SRASegment]:
        return [s for s in self.segments if s.read_type == "Biological"]

    def technical_segments(self) -> list[SRASegment]:
        return [s for s in self.segments if s.read_type == "Technical"]

    def get_segment(self, index: int) -> SRASegment:
        return self.segments[index]

    def get_sequence(self, index: int) -> str:
        return self.segments[index].sequence

    def get_quality(self, index: int) -> str:
        return self.segments[index].quality


class SRAReader:
    """
    Streaming reader for SRA files.

    On init, detects read structure (lengths, starts, types) from metadata.
    If structure is fixed (the common case for Illumina), uses a Cython
    fast path that only fetches the READ column and slices at C speed
    with precomputed offsets. No CSV parsing per line.

    Usage:
        reader = SRAReader("SRR12345678.sra")
        for spot in reader:
            for seg in spot.segments:
                process(seg.sequence)

        # Fastest — raw tuples of bytes, no Python objects:
        for segs in reader.iter_raw():
            barcode, cdna = segs[0], segs[1]
    """

    def __init__(
        self,
        sra_path: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        include_quality: bool = True,
        info: Optional[SRAInfo] = None,
    ):
        check_dependencies()
        self.sra_path = sra_path
        self.start = start
        self.end = end
        self.include_quality = include_quality

        self._info = info or sra_info(sra_path)

        # Precompute segment slicing from fixed structure
        self._seg_slices: tuple[tuple[int, int], ...] = ()
        self._seg_types: list[str] = []
        self._is_fixed = not self._info.is_variable_length
        self._precompute_structure()

        self._process: Optional[subprocess.Popen] = None

    def _precompute_structure(self):
        read_lengths = self._info.read_lengths
        read_types = self._info.read_types
        slices = []
        offset = 0
        for length, rtype in zip(read_lengths, read_types):
            slices.append((offset, offset + length))
            self._seg_types.append(rtype)
            offset += length
        self._seg_slices = tuple(slices)

    @property
    def info(self) -> SRAInfo:
        return self._info

    @property
    def num_spots(self) -> int:
        return self._info.num_spots

    def _row_range(self) -> Optional[tuple[int, int]]:
        if self.start is not None or self.end is not None:
            return (self.start or 1, self.end or self._info.num_spots)
        return None

    # ------------------------------------------------------------------
    # Iteration modes (fastest → richest)
    # ------------------------------------------------------------------

    def iter_raw(self, backend: Optional[str] = None) -> Iterator[tuple[bytes, ...]]:
        """
        Fastest iteration: yields tuples of segment byte strings.

        Auto-selects backend:
        - "xsra": Uses xsra (VDB C library, multi-threaded). ~5-20x faster.
        - "vdb_dump": Uses vdb-dump subprocess + Cython parsing.

        If backend is None, uses xsra if available, else vdb-dump.

        Yields:
            (seg0_bytes, seg1_bytes, ...) per spot.
        """
        if backend is None:
            backend = "xsra" if has_xsra() else "vdb_dump"

        if backend == "xsra":
            yield from self._iter_raw_xsra()
        else:
            yield from self._iter_raw_vdb_dump()

    def _iter_raw_vdb_dump(self) -> Iterator[tuple[bytes, ...]]:
        """vdb-dump backend for iter_raw."""
        proc = run_vdb_dump(self.sra_path, ["READ"], self._row_range(), fmt="csv")
        self._process = proc
        slices = self._seg_slices

        try:
            for raw_line in proc.stdout:
                line = raw_line.rstrip(b"\n")
                if not line:
                    continue
                yield parse_line_fixed_seqonly(line, slices)
        finally:
            proc.stdout.close()
            proc.wait()
            self._process = None

    def _iter_raw_xsra(self, n_threads: int = 4) -> Iterator[tuple[bytes, ...]]:
        """
        xsra backend: reads SRA via the VDB C library directly.

        xsra outputs FASTA to stdout. With --split-spot, segments for the
        same spot appear consecutively. The output is:
            >spot_id.seg_idx
            SEQUENCE
            >spot_id.seg_idx
            SEQUENCE
            ...

        We group consecutive segments by spot to yield tuples.
        """
        num_segs = len(self._seg_slices)
        seg_indices = list(range(num_segs))

        limit = None
        rr = self._row_range()
        if rr is not None:
            limit = rr[1] - rr[0] + 1

        proc = run_xsra_stream(
            self.sra_path,
            segment_indices=seg_indices,
            row_limit=limit,
            n_threads=n_threads,
            fmt="fasta",
        )
        self._process = proc

        try:
            segments: list[bytes] = []
            for raw_line in proc.stdout:
                line = raw_line.rstrip(b"\n")
                if not line:
                    continue
                if line.startswith(b">"):
                    # Header line — skip, but if we've accumulated enough
                    # segments for a full spot, yield it
                    pass
                else:
                    # Sequence line
                    segments.append(line)
                    if len(segments) == num_segs:
                        yield tuple(segments)
                        segments = []
            # Yield any remaining partial spot
            if segments:
                yield tuple(segments)
        finally:
            proc.stdout.close()
            proc.wait()
            self._process = None

    def iter_raw_batched(
        self, batch_size: int = 10000, backend: Optional[str] = None,
    ) -> Iterator[tuple[list[bytes], ...]]:
        """
        Batch iteration: yields M lists (one per segment) of N byte strings.

        Ideal for feeding into vectorized Cython kmer processing.
        Uses xsra when available (5-20x faster than vdb-dump).

        Args:
            batch_size: Number of spots per batch.
            backend: ``"xsra"`` or ``"vdb_dump"``. Auto-selects xsra when available.

        Yields:
            ([seg0_spot0, seg0_spot1, ...], [seg1_spot0, ...], ...)
        """
        if backend is None:
            backend = "xsra" if has_xsra() else "vdb_dump"

        if backend == "xsra":
            yield from self._iter_raw_batched_xsra(batch_size)
        else:
            yield from self._iter_raw_batched_vdb_dump(batch_size)

    def _iter_raw_batched_xsra(
        self, batch_size: int,
    ) -> Iterator[tuple[list[bytes], ...]]:
        """xsra backend for iter_raw_batched."""
        num_segs = len(self._seg_slices)
        batch: list[list[bytes]] = [[] for _ in range(num_segs)]
        count = 0

        for spot_segs in self._iter_raw_xsra():
            for i, seq in enumerate(spot_segs):
                batch[i].append(seq)
            count += 1
            if count >= batch_size:
                yield tuple(batch)
                batch = [[] for _ in range(num_segs)]
                count = 0

        if count > 0:
            yield tuple(batch)

    def _iter_raw_batched_vdb_dump(
        self, batch_size: int,
    ) -> Iterator[tuple[list[bytes], ...]]:
        """vdb-dump backend for iter_raw_batched."""
        from ._parser import parse_lines_batch_seqonly

        proc = run_vdb_dump(self.sra_path, ["READ"], self._row_range(), fmt="csv")
        self._process = proc
        slices = self._seg_slices

        try:
            batch: list[bytes] = []
            for raw_line in proc.stdout:
                line = raw_line.rstrip(b"\n")
                if not line:
                    continue
                batch.append(line)
                if len(batch) >= batch_size:
                    yield parse_lines_batch_seqonly(batch, slices)
                    batch = []
            if batch:
                yield parse_lines_batch_seqonly(batch, slices)
        finally:
            proc.stdout.close()
            proc.wait()
            self._process = None

    def __iter__(self) -> Iterator[SRASpot]:
        """
        Standard iteration: yields SRASpot objects with SRASegment members.

        Uses Cython fast path for fixed structure, Python fallback otherwise.
        """
        if self._is_fixed:
            yield from self._iter_fast()
        else:
            yield from self._iter_full()

    def _iter_fast(self) -> Iterator[SRASpot]:
        """Fast path: fixed structure, Cython line parsing."""
        proc = run_vdb_dump(
            self.sra_path,
            ["READ", "QUALITY"] if self.include_quality else ["READ"],
            self._row_range(),
            fmt="csv",
        )
        self._process = proc
        slices = self._seg_slices
        seg_types = self._seg_types
        has_qual = self.include_quality
        spot_id = self.start or 1

        try:
            for raw_line in proc.stdout:
                line = raw_line.rstrip(b"\n")
                if not line:
                    continue

                parsed = parse_line_fixed(line, slices, has_qual)

                segments = []
                for i, item in enumerate(parsed):
                    if has_qual:
                        seq_bytes, qual_bytes = item
                        seq = seq_bytes.decode("ascii")
                        qual = qual_bytes.decode("ascii")
                    else:
                        seq = item.decode("ascii")
                        qual = ""
                    segments.append(SRASegment(
                        sequence=seq,
                        quality=qual,
                        read_type=seg_types[i],
                        segment_index=i,
                    ))
                yield SRASpot(spot_id=spot_id, segments=segments)
                spot_id += 1
        finally:
            proc.stdout.close()
            proc.wait()
            self._process = None

    def _iter_full(self) -> Iterator[SRASpot]:
        """Slow fallback: variable-length reads, full CSV parsing per line."""
        columns = ["READ", "READ_TYPE", "READ_LEN", "READ_START"]
        if self.include_quality:
            columns.append("QUALITY")

        proc = run_vdb_dump(self.sra_path, columns, self._row_range(), fmt="csv")
        self._process = proc
        spot_id = self.start or 1

        try:
            for raw_line in proc.stdout:
                line = raw_line.rstrip(b"\n").decode("ascii", errors="replace")
                if not line:
                    continue
                try:
                    fields = _parse_vdb_csv_line(line)
                    if len(fields) < 4:
                        spot_id += 1
                        continue

                    full_seq = fields[0]
                    read_types = [_parse_read_type(t) for t in fields[1].split(",")]
                    read_lens = [int(l.strip()) for l in fields[2].split(",")]
                    read_starts = [int(s.strip()) for s in fields[3].split(",")]
                    full_qual = fields[4] if self.include_quality and len(fields) > 4 else ""

                    segments = []
                    for i, (rtype, rlen, rstart) in enumerate(
                        zip(read_types, read_lens, read_starts)
                    ):
                        segments.append(SRASegment(
                            sequence=full_seq[rstart:rstart + rlen],
                            quality=full_qual[rstart:rstart + rlen] if full_qual else "",
                            read_type=rtype,
                            segment_index=i,
                        ))
                    yield SRASpot(spot_id=spot_id, segments=segments)
                    spot_id += 1
                except (ValueError, IndexError):
                    spot_id += 1
                    continue
        finally:
            proc.stdout.close()
            proc.wait()
            self._process = None

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    def iter_segment(self, segment_index: int) -> Iterator[tuple[bytes, bytes]]:
        """Iterate yielding (sequence, quality) bytes for one segment."""
        slices = self._seg_slices
        s, e = slices[segment_index]

        proc = run_vdb_dump(
            self.sra_path,
            ["READ", "QUALITY"] if self.include_quality else ["READ"],
            self._row_range(),
            fmt="csv",
        )
        self._process = proc

        try:
            for raw_line in proc.stdout:
                line = raw_line.rstrip(b"\n")
                if not line:
                    continue
                if self.include_quality:
                    parsed = parse_line_fixed(line, slices, True)
                    yield parsed[segment_index]
                else:
                    yield (line[s:e], b"")
        finally:
            proc.stdout.close()
            proc.wait()
            self._process = None

    def iter_sequences_only(self) -> Iterator[tuple[bytes, ...]]:
        """Alias for iter_raw()."""
        return self.iter_raw()

    def iter_raw_parallel(
        self, n_workers: int = 4, queue_size: int = 10000,
    ) -> Iterator[tuple[bytes, ...]]:
        """
        Parallel iteration: spawns N vdb-dump processes on row-range chunks.

        Each worker reads a disjoint range of spots. A background thread per
        worker pushes parsed tuples into a shared queue. The main thread
        yields from the queue in chunk order (worker 0 fully, then worker 1, etc.)
        to preserve spot ordering.

        This saturates the CPU across multiple vdb-dump decompression threads,
        which is the actual bottleneck at ~96k spots/sec single-threaded.

        Args:
            n_workers: Number of parallel vdb-dump processes.
            queue_size: Max items buffered per worker before backpressure.

        Yields:
            (seg0_bytes, seg1_bytes, ...) per spot, in order.
        """
        import threading
        from queue import Queue

        total = self._info.num_spots
        s_start = self.start or 1
        s_end = self.end or total

        span = s_end - s_start + 1
        chunk_size = span // n_workers
        if chunk_size < 1:
            yield from self.iter_raw()
            return

        slices = self._seg_slices
        _sentinel = None

        def _worker(worker_range, out_q):
            """Read a range of spots and push parsed tuples to queue."""
            proc = run_vdb_dump(
                self.sra_path, ["READ"], worker_range, fmt="csv",
            )
            try:
                for raw_line in proc.stdout:
                    line = raw_line.rstrip(b"\n")
                    if line:
                        out_q.put(parse_line_fixed_seqonly(line, slices))
            finally:
                proc.stdout.close()
                proc.wait()
                out_q.put(_sentinel)

        # Build ranges and launch workers in order
        queues = []
        threads = []
        for i in range(n_workers):
            ws = s_start + i * chunk_size
            we = s_start + (i + 1) * chunk_size - 1 if i < n_workers - 1 else s_end
            q = Queue(maxsize=queue_size)
            queues.append(q)
            t = threading.Thread(target=_worker, args=((ws, we), q), daemon=True)
            threads.append(t)
            t.start()

        # Yield in order: drain worker 0 fully, then worker 1, etc.
        for q in queues:
            while True:
                item = q.get()
                if item is _sentinel:
                    break
                yield item

        for t in threads:
            t.join()

    def close(self):
        if self._process is not None:
            self._process.terminate()
            self._process.wait()
            self._process = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def chunk_ranges(self, num_chunks: int) -> list[tuple[int, int]]:
        """Divide spots into N ranges for parallel processing."""
        total = self._info.num_spots
        chunk_size = total // num_chunks
        ranges = []
        for i in range(num_chunks):
            s = i * chunk_size + 1
            e = (i + 1) * chunk_size if i < num_chunks - 1 else total
            ranges.append((s, e))
        return ranges
