# Mitogenome Coverage Pipeline

Per-base read coverage for circularized mitochondrial genomes, corrected for the artificial coverage dip at the point where the circle is cut to write a linear FASTA file.

**Version 4.1** (November 2026) · Bowtie2 · SAMtools · Python (NumPy, pandas, Matplotlib)

## Why a correction is needed

A circular mitogenome is a closed loop, but its FASTA file is linear, with an arbitrary cut point between the last base and the first. Aligners treat that cut point as a real end. Reads that span it cannot align contiguously, and read pairs with one mate on either side look a whole genome apart, so they fail proper-pair filtering. The result is a coverage dip at both ends of the sequence that is easily mistaken for a poorly assembled region.

The pipeline pads the reference by appending its first *K* bp to its 3′ end, where *K* is the mean insert size plus four standard deviations (minimum 300 bp), so that origin-spanning pairs have somewhere contiguous to map. It then folds the depth from the appended region back onto positions 1 to *K* before computing statistics. A control mapping to the unpadded assembly checks that the correction does not alter depth anywhere else.

## Requirements

| Item | Version or format | Install or notes |
| --- | --- | --- |
| Bowtie2 | 2.4.4 used in our analyses | `sudo apt install bowtie2` |
| SAMtools | 1.13 | `sudo apt install samtools` |
| Python | 3.8 or later | `pip install numpy pandas matplotlib` |
| Assembly | FASTA, one circularized contig | Nucleotide and IUPAC codes only (no NOVOPlasty `*` flags) |
| Reads | Paired-end FASTQ (R1, R2) | Single-end data are not supported |

The mapping steps run under Linux (on Windows, inside WSL with Ubuntu). The final Python step also runs in Windows Python/IDLE.

## Quick start

Set up one folder per specimen:

```
<specimen>/
    assembly.fasta            circularized mitogenome (rotated, if you rotated it)
    original_assembly.fasta   if rotated: the unrotated assembler output
    R1.fastq
    R2.fastq
    scripts/                  the contents of this repository's scripts/ folder
```

Open `scripts/run_pipeline.sh` and check the settings block at the top. Set `THREADS`, and if you rotated the assembly (for example to begin at *cox1*), set `ORIGINAL_ASSEMBLY="original_assembly.fasta"` so the driver verifies the rotation.

1. Check the assembler output for leftover circularization overlap, before any rotation (exit code 0 = none of 20 bp or more):
   ```
   python3 scripts/00_check_terminal_overlap.py original_assembly.fasta
   ```
2. Run the driver from the specimen folder. Everything lands in `coverage_run/`:
   ```
   bash scripts/run_pipeline.sh
   ```
3. Copy `padded_depth.txt`, `padding_metadata.txt` and `prelim_filtered_depth.txt` from `coverage_run/` next to `02_fold_and_plot.py`.
4. Set `PLOT_TITLE` near the top of `02_fold_and_plot.py`, then run it and enter an output prefix (for example `Genus_species_specimenID`) when asked:
   ```
   python3 02_fold_and_plot.py
   ```

## Outputs

| File | Contents | What to look for |
| --- | --- | --- |
| `<prefix>_coverage_stats.txt` | Genome length, *K*, mean, SD, minimum, maximum, breadth of coverage | The values to report |
| `<prefix>_control_comparison.txt` | Corrected vs unpadded depth across the genome interior (positions *K*+1 to *N*−*K*) | Difference 0% or very close; same minimum at the same position(s) |
| `<prefix>_linear_coverage.png`, `<prefix>_circular_coverage.png` | Coverage figures (300 dpi) | No dip at either end of the genome |
| `coverage_run/run_provenance.txt` | Date, software versions, settings, input checksums, check results | Keep with the results |

The interior means in the comparison file leave out *K* bp at each end, so they differ from the whole-genome mean; report the whole-genome values from the statistics file.

## Validation

Padding changes only how reads near the origin are placed, so alignments in the genome interior are unaffected and the corrected and unpadded profiles should match there. In every specimen tested so far, interior mean depth agreed exactly or within 0.5%, with the same minimum at the same positions, across mean depths from about 40× to about 30,700×.

## Provenance

Each run records how it was produced. `run_pipeline.sh` writes `run_provenance.txt` (software versions, settings, MD5 checksums of every input and script, check results, *K*, alignment counts). `padding_metadata.txt` records the assembly checksum, and the statistics and comparison files end with a provenance block giving the script version, date, Python library versions and input checksums. Mapping commands are also stored in each BAM header (`samtools view -H filtered.bam | grep "^@PG"`).

## Important notes

- **Do not add a MAPQ filter (`-q`).** Reads in the duplicated padding window align equally well to both copies and receive low mapping quality by design. Filtering them reintroduces the artifact. The proper-pair filter (`-f 2`) provides the quality control.
- **Use Bowtie2 in `--local` mode**, as the driver does. End-to-end mode is unforgiving at the padded junction.
- **Rotated assemblies:** run the overlap check on the assembler output before rotating, then verify that the rotated file is an exact cyclic shift (`verify_rotation.py`, or the `ORIGINAL_ASSEMBLY` setting).

## Repository contents

| Path | Purpose |
| --- | --- |
| `scripts/run_pipeline.sh` | Driver: checks, mapping, padding, filtering, depth, control depth, provenance record |
| `scripts/00_check_terminal_overlap.py` | Checks for leftover circularization overlap and invalid characters |
| `scripts/verify_rotation.py` | Confirms a rotated assembly is an exact cyclic shift of the original |
| `scripts/01_build_padded_reference.py` | Builds the padded reference and its metadata |
| `scripts/02_fold_and_plot.py` | Folds depth, computes statistics, draws figures, compares with the control |
| `scripts/03_compare_to_control.py` | Optional standalone control comparison, for profiles folded with earlier versions |
| `scripts/diagnose_padding_artifact.sh` | Optional troubleshooting report for the padded window |
| `docs/Coverage_pipeline_quick_start.docx` | Two-page quick start guide |
| `docs/Coverage_calculation_and_visualization_protocol_v4.1.docx` | Full step-by-step protocol with embedded scripts |
| `docs/README_full_reference.txt` | Full reference: settings, outputs, troubleshooting, limitations |
| `docs/Supplement_1_command_line_protocol.docx` | Condensed command-line protocol for manuscript supplements |
| `docs/Coverage_Methods.docx` | Methods and Results template for manuscripts |
| `CHANGELOG.md` | Version history |

## Limitations

The pipeline expects a single circular contig and paired-end reads. The overlap check finds only exact suffix–prefix matches, so an overlap containing mismatches is not detected, and a genuine repeat at the origin can trigger it. Folding assumes one primary alignment per read, which holds with Bowtie2's default reporting (not with `-k`).

## References

1. Dierckxsens N, Mardulyn P, Smits G. 2017. NOVOPlasty: de novo assembly of organelle genomes from whole genome data. *Nucleic Acids Research* 45(4): e18. [doi:10.1093/nar/gkw955](https://doi.org/10.1093/nar/gkw955)
2. Langmead B, Salzberg SL. 2012. Fast gapped-read alignment with Bowtie 2. *Nature Methods* 9(4): 357–359. [doi:10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)
3. Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. 2021. Twelve years of SAMtools and BCFtools. *GigaScience* 10(2): giab008. [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)
4. Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, Cournapeau D, Wieser E, Taylor J, Berg S, Smith NJ, Kern R, Picus M, Hoyer S, van Kerkwijk MH, Brett M, Haldane A, Fernández del Río J, Wiebe M, Peterson P, Gérard-Marchant P, Sheppard K, Reddy T, Weckesser W, Abbasi H, Gohlke C, Oliphant TE. 2020. Array programming with NumPy. *Nature* 585(7825): 357–362. [doi:10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)
5. Hunter JD. 2007. Matplotlib: a 2D graphics environment. *Computing in Science & Engineering* 9(3): 90–95. [doi:10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55)
6. McKinney W. 2010. Data structures for statistical computing in Python. In: van der Walt S, Millman J, editors. *Proceedings of the 9th Python in Science Conference*. p. 56–61. [doi:10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a)
