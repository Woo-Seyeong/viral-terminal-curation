# viral-terminal-curation
Python/Jupyter workflow for curating hairpin-ended dsDNA virus genome termini — DRC trimming, orientation normalization, and ITR candidate detection.


## Overview

Hairpin-ended dsDNA viruses such as chloroviruses carry structurally complex genome termini (hairpin ends and inverted terminal repeats, ITRs) that are often misrepresented in genome assemblies. A common artifact is **duplicated reverse-complementary (DRC) sequence**, in which redundant terminal sequence is included in reverse-complement orientation at a genome end. If left uncorrected, DRC artifacts can inflate genome length, duplicate or fragment terminal genes, and interfere with ITR detection.

This repository provides a reproducible workflow that detects and removes DRC artifacts, refines terminal boundaries, normalizes genome orientation against a reference, and reports ITR candidates, while preserving the evidence behind each automated decision for manual review.

The workflow was developed as part of an M.S. thesis on terminal-structure analysis of *Chlorovirus vanettense* PKV2 and 130 public chlorovirus genome records.

## What the workflow does

Each input sequence passes through four main steps:

1. **DRC detection and trimming** — a coarse scan locates a candidate terminal region, and a fine-scale transition-run check (opposite → overlap_or_none → reversed states) validates a true DRC boundary before trimming.
2. **Terminal boundary refinement** — residual hairpin- or palindrome-associated sequence remaining after primary DRC trimming is removed.
3. **Reference-guided orientation normalization** — genome orientation is normalized using k-mer support against a reference genome (e.g., PBCV-1), when a reference is provided.
4. **ITR candidate detection** — left and reverse-complemented right termini are compared by local alignment and filtered to report ITR candidates.

## Requirements

- Python 3.10.20
- Biopython 1.87
- Jupyter Notebook
- (standard scientific Python stack)

```bash
pip install biopython notebook
```

## Usage

1. Open the notebook:
```bash
   jupyter notebook hairpin_finder_v2.ipynb
```
2. Edit the **`# User configuration`** cell at the top to set input/output paths, the reference genome (for orientation normalization), and parameters.
3. Run all cells.

The workflow accepts genome sequences in FASTA format and supports both single-record and multi-record FASTA files.

## Parameters

All parameters are defined in the configuration cell at the top of the notebook. Defaults used in the thesis include:

| Parameter | Value |
|---|---|
| DRC window size | 300 bp |
| DRC coarse / fine step | 2,000 bp / 100 bp |
| DRC min identity / coverage | 0.85 / 0.80 |
| DRC margin | 300 bp |
| Max ambiguous-base fraction | 0.20 |
| Transition run (opposite / reversed) | ≥2 / ≥2 |
| Orientation k-mer size | 15 |
| Orientation support ratio | 1.1 |
| Min ITR length / identity | 100 bp / 0.95 |
| ITR terminal tolerance | 500 bp |

## Output files

| File | Contents |
|---|---|
| `*_curated.fasta` | Curated sequence for each input genome |
| combined FASTA | All curated sequences, for downstream comparative analysis |
| `summary_report.csv` | Genome-level curation results (lengths, DRC cuts, orientation status, ITR candidate info) |
| `drc_debug_report.csv` | Window-level evidence supporting each DRC trimming decision |

## Citation

If you use this workflow, please cite the associated thesis (see `CITATION.cff`).

## License

Released under the MIT License. See `LICENSE`.
