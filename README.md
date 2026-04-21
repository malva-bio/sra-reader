# sra-reader

Direct SRA file reader for single-cell experiments. Bypasses `fasterq-dump` by streaming reads directly from `.sra` files with a Cython-optimized parsing loop.

## Install

```bash
pip install .
```

Requires `sra-tools` on PATH:
```bash
conda install -c bioconda sra-tools
```

## Usage

```python
from sra_reader import SRAReader, detect_layout, sra_info

# Inspect structure
info = sra_info("SRR12345678.sra")
layout = detect_layout("SRR12345678.sra", info=info)
print(layout.describe())

# Fastest: raw bytes tuples (for Cython kmer indexers)
reader = SRAReader("SRR12345678.sra", include_quality=False, info=info)
for segs in reader.iter_raw():
    barcode, cdna = segs[0], segs[1]

# Batched: per-segment lists (for vectorized processing)
for seg_lists in reader.iter_raw_batched(batch_size=50000):
    seg0_batch, seg1_batch = seg_lists  # list[bytes] each

# Rich objects (when you need metadata)
for spot in reader:
    for seg in spot.segments:
        print(seg.read_type, seg.sequence[:20])
```
