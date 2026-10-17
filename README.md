# viral-terminal-curation

Python/Jupyter workflow for curating hairpin-ended dsDNA virus genome termini — DRC trimming, orientation normalization, and ITR candidate detection.

## Open in Colab

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Woo-Seyeong/viral-terminal-curation/blob/main/hairpin_terminal_curator_v3.1.ipynb)

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

## Repository contents

| File | Purpose |
|---|---|
| `hairpin_terminal_curator_v3.1.ipynb` | The curation workflow: DRC detection and trimming, terminal boundary refinement, orientation normalization, ITR candidate detection. |
| `simulation_benchmark.py` | Generates the simulated hairpin-ended genomes and reproduces the benchmark of the accompanying paper (Table 3, Figure S1). |
| `PKV2_uncurated.fasta` | The uncurated *Chlorovirus vanettense* PKV2 long-read assembly (398,401 bp), the input from which the curated genome was derived. |
| `illumina_contigs/` | Short-read (Illumina) assembly contigs for PKV2 used to check the long-read consensus. |

## Data

### `PKV2_uncurated.fasta`

Single Canu v1.8 contig (`tig00000001`, 398,401 bp, 17,453 reads) assembled from Oxford Nanopore PromethION reads of PKV2. Both termini carry duplicated reverse-complementary (DRC) sequence. This file is the input used for every PKV2 result in the accompanying paper: the manual terminal curation, the automated workflow run, and the assembler-agreement comparison. Running the notebook on it reproduces the curated 309,183-bp genome deposited as GenBank **PZ278399.1**.

The deposited record is the curated sequence, so this uncurated assembly is not retrievable from GenBank and is provided here instead.

### `illumina_contigs/`

Contigs from two independent short-read assemblies of PKV2, used in the original analysis to check the long-read consensus base by base. The prefixes `20-30_` and `30-40_` are the sample labels used in that study and carry no other meaning.

Both samples were prepared as TruSeq Nano DNA libraries and sequenced as 2 × 151-bp paired-end reads on an Illumina NovaSeq 6000 (Macrogen, Seoul, Republic of Korea, March 2022). Reads were trimmed with Trimmomatic v0.36, randomly subsampled to 298,014 reads (≈45 Mb, 146× coverage) and assembled with Platanus-allee v2.2.2. The raw reads were held by the sequencing provider for three months only and are no longer available.

Aligned against the curated genome (`minimap2 -x asm5`), the two sets agree with the deposited sequence to within a handful of substituted positions; the large gaps are collapsed repeats, as expected of a short-read assembly:

| Set | Contigs | Assembly size (N50) | Aligned range on PZ278399.1 | Substitutions | Left end not covered | Right end not covered |
|---|---|---|---|---|---|---|
| `20-30` | 4 | 305,602 bp (229,885) | 25–306,616 | 3 | 24 bp | 2,567 bp |
| `30-40` | 3 | 305,394 bp (242,680) | 2,366–309,028 | 12 | 2,365 bp | 155 bp |

Neither set reaches either genome end, so the terminal hairpin loops are not fully covered by short-read data. One contig does, however, extend into the left-terminal loop: `20-30_contig2.fa` begins at position 25 of the curated genome and its only two mismatches against the deposited sequence over positions 25–69 are at **positions 32 and 35**. Substituting those two bases reproduces exactly the 37-bp left loop sequence carried by 95% of the nanopore reads that cross it:

```
PZ278399.1  (as deposited)  TAGCTATACAGCCTGTTTAAAGTGGCAATTTAGCTACG
nanopore reads / contig     TAGCTATACAGCCTGTTTAAAGTGGCAATTTCGCGAC
                                                           ^  ^
                                                          32 35
```

Both positions lie within the terminal 11 bp of that contig, which a global aligner soft-clips, which is why the discrepancy was not caught in the original alignment-based correction pass. The deposited record is being corrected at these two positions; the correction changes no other base, and neither genome length nor any trimming boundary depends on it.

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
   jupyter notebook hairpin_terminal_curator_v3.1.ipynb
   ```
2. In the configuration cell, set `INPUT_FASTA` and `OUTPUT_DIR` (and optionally `REFERENCE_FASTA`).
3. Run all cells.

### Google Colab

1. Click the **Open in Colab** badge above, and run the first cell to install Biopython.
2. Upload your FASTA using the left **Files** panel (drag & drop), then set `INPUT_FASTA` in the configuration cell to its path (e.g. `/content/Input_file.fasta`).
3. **Runtime ▸ Run all**.

## Parameters

All parameters are defined in the configuration cell at the top of the notebook. Defaults used in the accompanying paper include:

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

## Versions

Releases are tagged. Results in the accompanying paper were produced with the workflow as tagged **v3.1** (commit `764af03`); cite that tag, or the commit it points to, rather than `main`. **v3.2** adds the data files above without changing the workflow.

## Citation

If you use this workflow, please cite the accompanying paper (details to be added on publication) and this repository.
