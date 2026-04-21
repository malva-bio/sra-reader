# sra-reader

Direct SRA file reader for single-cell experiments. Bypasses `fasterq-dump` by streaming reads directly from `.sra` files with a Cython-optimized parsing loop.

## Install

```bash
pip install .
```

### Required: sra-tools

`vdb-dump` and `sra-stat` must be on your PATH — they do the actual SRA decoding:

```bash
conda install -c bioconda sra-tools
```

### Optional: xsra (recommended for best performance)

`xsra` is an alternative SRA backend from the Arc Institute that reads via the VDB C library directly, bypassing text serialization. when available, `iter_raw()` picks it up automatically and runs **5–20× faster** than `vdb-dump`. it is not required — the reader falls back to `vdb-dump` if `xsra` is not on PATH.

install via cargo (requires a Rust toolchain):

```bash
cargo install xsra
```

if you don't have Rust:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cargo install xsra
```

> **note**: only `iter_raw()` uses `xsra`. `iter_raw_batched()` always uses `vdb-dump`.

## Usage

```python
from sra_reader import SRAReader, detect_layout, sra_info

# inspect structure
info = sra_info("SRR12345678.sra")
layout = detect_layout("SRR12345678.sra", info=info)
print(layout.describe())

# fastest: raw bytes tuples — uses xsra if available, vdb-dump otherwise
reader = SRAReader("SRR12345678.sra", include_quality=False, info=info)
for segs in reader.iter_raw():
    barcode, cdna = segs[0], segs[1]

# force a specific backend
for segs in reader.iter_raw(backend="xsra"):    # requires xsra on PATH
    ...
for segs in reader.iter_raw(backend="vdb_dump"):  # always available
    ...

# batched: per-segment lists, good for vectorized processing (always vdb-dump)
for seg_lists in reader.iter_raw_batched(batch_size=50000):
    seg0_batch, seg1_batch = seg_lists  # list[bytes] each

# rich objects: when you need read types, quality, etc.
for spot in reader:
    for seg in spot.segments:
        print(seg.read_type, seg.sequence[:20])
```
