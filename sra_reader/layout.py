"""
Read layout detection for single-cell experiments.

Samples N spots and uses heuristics on read lengths and sequence
composition to figure out which segment is barcode vs cDNA.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
from collections import Counter

from .reader import SRAReader, SRASpot
from .utils import SRAInfo, sra_info


KNOWN_BARCODE_LENGTHS = {
    28: "10x_v2_or_v3",
    26: "10x_v2",
    20: "dropseq",
    24: "indrops",
    16: "CEL-seq2",
    9:  "sci-RNA-seq",
    10: "sci-RNA-seq",
}

KNOWN_INDEX_LENGTHS = {6, 8, 10}


@dataclass
class ReadLayout:
    """Describes the read layout of an SRA accession."""
    barcode_segment: Optional[int]
    cdna_segment: Optional[int]
    index_segments: list[int]
    umi_segment: Optional[int]
    barcode_length: Optional[int]
    umi_length: Optional[int]
    umi_start: Optional[int]
    platform_guess: str
    num_segments: int
    segment_labels: list[str]
    confidence: float

    @property
    def is_paired_end(self) -> bool:
        return (
            self.barcode_segment is not None
            and self.cdna_segment is not None
        )

    def describe(self) -> str:
        lines = [
            f"Platform guess: {self.platform_guess} (confidence: {self.confidence:.1%})",
            f"Segments per spot: {self.num_segments}",
        ]
        for i, label in enumerate(self.segment_labels):
            lines.append(f"  Segment {i}: {label}")
        if self.barcode_segment is not None:
            lines.append(
                f"Barcode: segment {self.barcode_segment}, "
                f"length {self.barcode_length}"
            )
        if self.umi_length:
            lines.append(
                f"UMI: segment {self.umi_segment}, "
                f"start {self.umi_start}, length {self.umi_length}"
            )
        if self.cdna_segment is not None:
            lines.append(f"cDNA: segment {self.cdna_segment}")
        return "\n".join(lines)


def detect_layout(
    sra_path: str,
    num_sample: int = 1000,
    info: Optional[SRAInfo] = None,
) -> ReadLayout:
    """
    Auto-detect the read layout of a single-cell SRA file.

    Samples spots and uses heuristics on read lengths, read type flags,
    and sequence composition to assign segment roles.
    """
    _info = info or sra_info(sra_path)
    spots = _sample_spots(sra_path, n=num_sample, info=_info)

    if not spots:
        return _unknown_layout(_info)

    num_segments = spots[0].num_segments

    seg_stats = [_SegmentStats() for _ in range(num_segments)]
    for spot in spots:
        for seg in spot.segments:
            i = seg.segment_index
            seg_stats[i].lengths.append(len(seg.sequence))
            seg_stats[i].types.append(seg.read_type)
            if seg.sequence:
                seg_stats[i].poly_t_fracs.append(
                    seg.sequence.count("T") / len(seg.sequence)
                )

    seg_summaries = [_summarize_segment(i, ss) for i, ss in enumerate(seg_stats)]
    return _assign_roles(seg_summaries, _info)


def _sample_spots(
    sra_path: str, n: int, info: SRAInfo,
) -> list[SRASpot]:
    """Read first N spots for layout detection."""
    n = min(n, info.num_spots)
    reader = SRAReader(
        sra_path, start=1, end=n, include_quality=False, info=info,
    )
    spots = []
    with reader:
        for spot in reader:
            spots.append(spot)
    return spots


@dataclass
class _SegmentStats:
    lengths: list[int] = field(default_factory=list)
    types: list[str] = field(default_factory=list)
    poly_t_fracs: list[float] = field(default_factory=list)


@dataclass
class _SegmentSummary:
    index: int
    median_length: int
    length_std: float
    is_fixed_length: bool
    dominant_type: str
    median_poly_t: float
    is_likely_index: bool
    is_likely_barcode: bool
    is_likely_cdna: bool


def _summarize_segment(index: int, stats: _SegmentStats) -> _SegmentSummary:
    lengths = stats.lengths
    if not lengths:
        return _SegmentSummary(
            index=index, median_length=0, length_std=0,
            is_fixed_length=True, dominant_type="Unknown",
            median_poly_t=0, is_likely_index=False,
            is_likely_barcode=False, is_likely_cdna=False,
        )

    sorted_lens = sorted(lengths)
    median_len = sorted_lens[len(sorted_lens) // 2]
    mean_len = sum(lengths) / len(lengths)
    variance = sum((l - mean_len) ** 2 for l in lengths) / max(len(lengths), 1)
    std = variance ** 0.5
    is_fixed = std < 2.0

    type_counts = Counter(stats.types)
    dominant = type_counts.most_common(1)[0][0] if type_counts else "Unknown"

    poly_t_sorted = sorted(stats.poly_t_fracs) if stats.poly_t_fracs else [0]
    median_poly_t = poly_t_sorted[len(poly_t_sorted) // 2]

    is_likely_index = (
        is_fixed and median_len in KNOWN_INDEX_LENGTHS and dominant == "Technical"
    )
    is_likely_barcode = (
        is_fixed and median_len in KNOWN_BARCODE_LENGTHS and not is_likely_index
    )
    is_likely_cdna = (
        (not is_fixed or median_len > 50)
        and not is_likely_index
        and not is_likely_barcode
    )

    return _SegmentSummary(
        index=index, median_length=median_len, length_std=std,
        is_fixed_length=is_fixed, dominant_type=dominant,
        median_poly_t=median_poly_t, is_likely_index=is_likely_index,
        is_likely_barcode=is_likely_barcode, is_likely_cdna=is_likely_cdna,
    )


def _assign_roles(
    summaries: list[_SegmentSummary], info: SRAInfo,
) -> ReadLayout:
    num_segments = len(summaries)
    labels = ["unknown"] * num_segments
    barcode_seg = None
    cdna_seg = None
    umi_seg = None
    index_segs = []

    for s in summaries:
        if s.is_likely_index:
            labels[s.index] = "index"
            index_segs.append(s.index)

    candidates = [s for s in summaries if labels[s.index] == "unknown"]
    barcode_candidates = [s for s in candidates if s.is_likely_barcode]
    if barcode_candidates:
        best = barcode_candidates[0]
        barcode_seg = best.index
        umi_seg = best.index
        labels[best.index] = "barcode+UMI"

    remaining = [s for s in candidates if labels[s.index] == "unknown"]
    cdna_candidates = [s for s in remaining if s.is_likely_cdna]
    if cdna_candidates:
        cdna_candidates.sort(key=lambda s: s.median_length, reverse=True)
        cdna_seg = cdna_candidates[0].index
        labels[cdna_seg] = "cDNA"

    remaining = [s for s in summaries if labels[s.index] == "unknown"]
    if barcode_seg is None and remaining:
        fixed_bio = [
            s for s in remaining
            if s.is_fixed_length and s.dominant_type == "Biological"
        ]
        if fixed_bio:
            fixed_bio.sort(key=lambda s: s.median_length)
            barcode_seg = fixed_bio[0].index
            umi_seg = fixed_bio[0].index
            labels[barcode_seg] = "barcode+UMI (inferred)"

    if cdna_seg is None:
        remaining = [s for s in summaries if labels[s.index] == "unknown"]
        if remaining:
            remaining.sort(key=lambda s: s.median_length, reverse=True)
            cdna_seg = remaining[0].index
            labels[cdna_seg] = "cDNA (inferred)"

    for s in summaries:
        if labels[s.index] == "unknown":
            labels[s.index] = f"unknown (len={s.median_length})"

    platform_guess, bc_len, umi_len, umi_start = _guess_platform(
        summaries[barcode_seg].median_length if barcode_seg is not None else 0
    )

    confidence = _compute_confidence(summaries, barcode_seg, cdna_seg)

    return ReadLayout(
        barcode_segment=barcode_seg, cdna_segment=cdna_seg,
        index_segments=index_segs, umi_segment=umi_seg,
        barcode_length=bc_len, umi_length=umi_len, umi_start=umi_start,
        platform_guess=platform_guess, num_segments=num_segments,
        segment_labels=labels, confidence=confidence,
    )


def _guess_platform(barcode_seg_length: int):
    if barcode_seg_length == 28:
        return "10x_v3", 16, 12, 16
    elif barcode_seg_length == 26:
        return "10x_v2", 16, 10, 16
    elif barcode_seg_length == 20:
        return "dropseq", 12, 8, 12
    elif barcode_seg_length == 24:
        return "indrops", 8, 6, 8
    elif barcode_seg_length in (9, 10):
        return "sci-RNA-seq", barcode_seg_length, None, 0
    elif barcode_seg_length == 16:
        return "CEL-seq2", 6, 6, 6
    elif barcode_seg_length == 0:
        return "unknown", None, None, None
    else:
        return "unknown_sc", barcode_seg_length, None, None


def _compute_confidence(summaries, barcode_seg, cdna_seg):
    score = 0.0
    total = 3.0
    if barcode_seg is not None:
        s = summaries[barcode_seg]
        if s.is_fixed_length and s.median_length in KNOWN_BARCODE_LENGTHS:
            score += 1.0
        elif s.is_fixed_length:
            score += 0.5
    if cdna_seg is not None:
        s = summaries[cdna_seg]
        if s.median_length > 50:
            score += 1.0
        elif s.median_length > 20:
            score += 0.5
    if barcode_seg is not None and cdna_seg is not None and barcode_seg != cdna_seg:
        score += 1.0
    return min(score / total, 1.0)


def _unknown_layout(info):
    return ReadLayout(
        barcode_segment=None, cdna_segment=None, index_segments=[],
        umi_segment=None, barcode_length=None, umi_length=None,
        umi_start=None, platform_guess="unknown",
        num_segments=info.num_reads_per_spot,
        segment_labels=["unknown"] * info.num_reads_per_spot,
        confidence=0.0,
    )
