# viral-terminal-curation

Python/Jupyter workflow for curating hairpin-ended dsDNA virus genome termini — DRC trimming, orientation normalization, and ITR candidate detection.

## Open in Colab

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Woo-Seyeong/viral-terminal-curation/blob/main/hairpin_finder_v2.ipynb)

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

- Python 3.10+
- biopython >= 1.80
- Jupyter Notebook (local use only)

```bash
pip install biopython notebook
```

On Google Colab no local installation is needed — the first notebook cell installs Biopython for you.

## Usage

Input/output paths are set in the **`# User configuration`** cell at the top of the notebook. The workflow accepts FASTA format and supports both single- and multi-record files. Set `REFERENCE_FASTA = None` to skip orientation normalization.

### Local (Jupyter)

1. Open the notebook:
   ```bash
   jupyter notebook hairpin_finder_v2.ipynb
   ```
2. In the configuration cell, set `INPUT_FASTA` and `OUTPUT_DIR` (and optionally `REFERENCE_FASTA`).
3. Run all cells.

### Google Colab

1. Click the **Open in Colab** badge above, and run the first cell to install Biopython.
2. Upload your FASTA using the left **Files** panel (drag & drop), then set `INPUT_FASTA` in the configuration cell to its path (e.g. `/content/Input_file.fasta`).
3. **Runtime ▸ Run all**.

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

All outputs are written under `OUTPUT_DIR`.

| File | Contents |
|---|---|
| `<header>/<header>.fasta` | Curated sequence for each input genome (one subfolder per genome) |
| `all_final_sequences.fasta` | All curated sequences, for downstream comparative analysis |
| `summary_report.csv` | Genome-level curation results (lengths, DRC cuts, orientation status, ITR candidate info) |
| `drc_debug_report.csv` | Window-level evidence supporting each DRC trimming decision |

## License

Released under the MIT License. See `LICENSE`.
