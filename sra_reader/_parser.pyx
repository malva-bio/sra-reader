# cython: language_level=3, boundscheck=False, wraparound=False
# cython: cdivision=True, initializedcheck=False
"""
Fast line parser for vdb-dump CSV output.

This is the hot path — called once per spot (millions of times).
Everything here operates on raw bytes to avoid Python object overhead.
"""

from libc.string cimport memchr, memcpy
from cpython.bytes cimport PyBytes_AS_STRING, PyBytes_GET_SIZE
from cpython.mem cimport PyMem_Malloc, PyMem_Free


def parse_line_fixed(
    bytes line,
    tuple seg_slices,
    bint has_quality,
):
    """
    Parse a vdb-dump line when read structure is fixed (the common case).

    When fetching READ only, each line is just the concatenated sequence:
        ACGTACGT...ACGTACGT

    When fetching READ,QUALITY, each line is:
        ACGTACGT...ACGTACGT,IIIIIIII...IIIIIIII

    Neither field ever contains commas, so finding the first comma is safe.

    Args:
        line: Raw bytes line from vdb-dump stdout (already stripped of newline).
        seg_slices: Tuple of (start, end) int pairs, precomputed from metadata.
        has_quality: Whether QUALITY column is present after a comma.

    Returns:
        Tuple of segment byte strings: (seg0_seq, seg1_seq, ...)
        Or if has_quality: ((seg0_seq, seg0_qual), (seg1_seq, seg1_qual), ...)
    """
    cdef:
        const char* buf = PyBytes_AS_STRING(line)
        Py_ssize_t buf_len = PyBytes_GET_SIZE(line)
        const char* seq_buf
        Py_ssize_t seq_len
        const char* qual_buf
        Py_ssize_t qual_len
        const char* comma_ptr
        Py_ssize_t i, s, e
        Py_ssize_t n_segs = len(seg_slices)

    if has_quality:
        # Find the comma separating READ from QUALITY
        comma_ptr = <const char*>memchr(buf, <int>b',', buf_len)
        if comma_ptr == NULL:
            # No comma found — treat entire line as sequence, no quality
            seq_buf = buf
            seq_len = buf_len
            qual_buf = NULL
            qual_len = 0
        else:
            seq_len = comma_ptr - buf
            seq_buf = buf
            qual_buf = comma_ptr + 1
            qual_len = buf_len - seq_len - 1
    else:
        seq_buf = buf
        seq_len = buf_len
        qual_buf = NULL
        qual_len = 0

    # Build output tuple — one entry per segment
    result = []
    for i in range(n_segs):
        s, e = seg_slices[i]

        # Bounds check the slice against actual sequence length
        if e > seq_len:
            e = seq_len
        if s > seq_len:
            s = seq_len

        seg_seq = line[s:e]

        if qual_buf != NULL and qual_len > 0:
            # Slice quality at same offsets
            qs = s
            qe = e
            if qe > qual_len:
                qe = qual_len
            if qs > qual_len:
                qs = qual_len
            # +1 offset because qual_buf starts after the comma
            seg_qual = line[(seq_len + 1 + qs):(seq_len + 1 + qe)]
            result.append((seg_seq, seg_qual))
        else:
            result.append(seg_seq)

    return tuple(result)


def parse_line_fixed_seqonly(
    bytes line,
    tuple seg_slices,
):
    """
    Fastest path: parse READ-only line into segment byte strings.

    No quality, no CSV, just slice a byte string at known offsets.
    This is the function your Cython kmer indexer should call.

    Args:
        line: Raw bytes line (just the concatenated sequence).
        seg_slices: Tuple of (start, end) int pairs.

    Returns:
        Tuple of bytes: (seg0_seq, seg1_seq, ...)
    """
    cdef:
        Py_ssize_t n_segs = len(seg_slices)
        Py_ssize_t buf_len = PyBytes_GET_SIZE(line)
        Py_ssize_t i, s, e

    result = []
    for i in range(n_segs):
        s, e = seg_slices[i]
        if e > buf_len:
            e = buf_len
        if s > buf_len:
            s = buf_len
        result.append(line[s:e])

    return tuple(result)


def parse_lines_batch_seqonly(
    list lines,
    tuple seg_slices,
):
    """
    Batch parse multiple lines into per-segment lists.

    Given N lines and M segments, returns M lists each of length N.
    This is ideal for feeding into vectorized Cython kmer processing.

    Args:
        lines: List of bytes lines.
        seg_slices: Tuple of (start, end) int pairs.

    Returns:
        Tuple of M lists: ([seg0_line0, seg0_line1, ...], [seg1_line0, ...], ...)
    """
    cdef:
        Py_ssize_t n_lines = len(lines)
        Py_ssize_t n_segs = len(seg_slices)
        Py_ssize_t i, j, s, e, buf_len
        bytes line

    # Pre-allocate output lists
    output = [[] for _ in range(n_segs)]
    for j in range(n_segs):
        (<list>output[j]).clear()

    for i in range(n_lines):
        line = <bytes>lines[i]
        buf_len = PyBytes_GET_SIZE(line)
        for j in range(n_segs):
            s, e = seg_slices[j]
            if e > buf_len:
                e = buf_len
            if s > buf_len:
                s = buf_len
            (<list>output[j]).append(line[s:e])

    return tuple(output)
