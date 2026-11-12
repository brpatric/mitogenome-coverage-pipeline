# Mitogenome Coverage Pipeline

Per-base read coverage for circularized mitochondrial genomes, corrected for the artificial coverage dip at the point where the circle is cut to write a linear FASTA file.

**Version 4.5** (November 2026) · Bowtie2 · SAMtools · Python (NumPy, pandas, Matplotlib)

## Why a correction is needed

A circular mitogenome is a closed loop, but its FASTA file is linear, with an arbitrary cut point between the last base and the first. Aligners treat that cut point as a real end. Reads that span it cannot align contiguously, and read pairs with one mate on either side look a whole genome apart, so they fail proper-pair filtering. The result is a coverage dip at both ends of the sequence that is easily mistaken for a poorly assembled region.

The pipeline pads the reference by appending its first *K* bp to its 3′ end, where *K* is the mean insert size plus four standard deviations (minimum 300 bp), so that origin-spanning pairs have somewhere contiguous to map. It then folds the depth from the appended region back onto positions 1 to *K* before computing statistics. Two checks follow. A control mapping to the unpadded assembly confirms that the correction does not alter depth in the genome interior. A second mapping, to a copy of the assembly cut at the opposite side of the circle, checks the ends of the sequence, where padding acts: if the folded profile departs from it there, the ends of the reported profile are taken from that mapping instead (see [End-window check](#end-window-check)).

## Requirements

| Item | Version or format | Install or notes |
| --- | --- | --- |
| Bowtie2 | 2.4.4 used in our analyses | `sudo apt install bowtie2` |
| SAMtools | 1.13 | `sudo apt install samtools` |
| Python | 3.8 or later | `pip install numpy pandas matplotlib` |
| Assembly | FASTA, one circularized contig from any organelle assembler (e.g. NOVOPlasty, Dierckxsens et al. 2017; GetOrganelle, Jin et al. 2020) | Nucleotide and IUPAC codes only; resolve any assembler-specific flags first (e.g. NOVOPlasty's `*`). If the assembler reports several sequences or paths, supply only the circular one you are using |
| Reads | Paired-end FASTQ (R1, R2), plain (`.fastq`) or gzipped (`.fastq.gz`) | Used where they are; no need to unzip or copy. Single-end data are not supported |

The mapping steps run under Linux (on Windows, inside WSL with Ubuntu). The final Python step also runs in Windows Python/IDLE.

## Quick start

Set up one folder per specimen:

```
<specimen>/
    assembly.fasta            circularized mitogenome (rotated, if you rotated it)
    original_assembly.fasta   if rotated: the unrotated assembler output
    R1.fastq.gz               forward reads (.fastq or .fastq.gz)
    R2.fastq.gz               reverse reads
    scripts/                  run_pipeline.sh, 00_check_terminal_overlap.py,
                              01_build_padded_reference.py, and (optional)
                              02_fold_and_plot.py
```

The file names above are the defaults; your files can have any names. Open `scripts/run_pipeline.sh` and edit the settings block at the top to match them:

```
ASSEMBLY="assembly.fasta"      # your mitogenome FASTA
R1="R1.fastq.gz"               # your forward reads
R2="R2.fastq.gz"               # your reverse reads
THREADS=8                      # processor threads available
ROTATED=""                     # "yes" or "no"; leave empty to be asked
ORIGINAL_ASSEMBLY=""           # if rotated: the unrotated assembler output
PREFIX=""                      # output prefix; empty = specimen folder name
PLOT_TITLE="Per-base Read Coverage"   # title on the coverage figures
```

- **`ASSEMBLY`, `R1`, `R2`:** set these to your mitogenome FASTA and your forward and reverse read files. Give each as a file name in the specimen folder, a path relative to it, or an absolute path, for example `Sample1_R1.fastq.gz`, `../shared_reads/Sample1_R1.fastq.gz` or `/mnt/c/data/Sample1_R1.fastq.gz`. Reads can be gzipped or not. The driver stops with a clear message if a file is not found.
- **`THREADS`:** the number of processor threads to use.
- **`ROTATED`, `ORIGINAL_ASSEMBLY`:** whether you rotated the assembly (for example to begin at *cox1*) and, if so, the unrotated assembler output, named the same way. If you leave them empty, the driver asks when it starts. A run with no one to answer, such as a cluster job, stops with a message instead of guessing, so set them there. For a rotated assembly the driver checks the unrotated file for leftover circularization overlap, verifies that the rotation is exact, and stops if either check fails.
- **`PREFIX`, `PLOT_TITLE`:** used only when `02_fold_and_plot.py` is in `scripts/` (see below). `PREFIX` names the output files; left empty, it is the specimen folder's name.
- Leave `MIN_PAD`, `SD_MULTIPLIER` and `MIN_OVERLAP` at their defaults (300, 4, 20).

**All in one run.** With `02_fold_and_plot.py` in `scripts/`, run the driver from the specimen folder and answer its questions if it asks:
```
bash scripts/run_pipeline.sh
```
The mapping files land in `coverage_run/`, and the statistics, figures, check results and report in `coverage_results/` (Phase 9). That is the whole procedure. The driver needs NumPy, pandas and Matplotlib installed in WSL for this step; if they are missing, it finishes the mapping, says so, and you can run `02_fold_and_plot.py` separately.

**In two parts.** Without `02_fold_and_plot.py` in `scripts/`, the driver stops after the mapping. Then:

1. Put `02_fold_and_plot.py` in any folder (WSL or Windows), and copy `padded_depth.txt`, `padding_metadata.txt`, `prelim_filtered_depth.txt` and `origin_check_depth.txt` from `coverage_run/` next to it. Copy `run_provenance.txt` too if you want the mapping details (software versions, insert size, rotation and overlap checks) in the report.
2. Set `PLOT_TITLE` near the top of `02_fold_and_plot.py`, then run it (F5 in IDLE, or from the shell) and enter an output prefix (for example `Genus_species_specimenID`) when asked:
   ```
   python3 02_fold_and_plot.py
   ```

### Optional steps

These are not needed for a routine run; `02_fold_and_plot.py` already runs both checks.

- **Standalone control comparison (`03_compare_to_control.py`).** Use this for a depth profile folded with a version of `02_fold_and_plot.py` earlier than 4.0, which did not run the comparison itself. Copy `<prefix>_folded_depth.txt`, `prelim_filtered_depth.txt` and `padding_metadata.txt` next to the script, then run it. It writes `<prefix>_control_comparison.txt` with the same values `02_fold_and_plot.py` would report:
  ```
  python3 03_compare_to_control.py
  ```
- **Padding diagnostic report.** The driver writes this automatically (Phase 6) as `coverage_run/diagnostic_report_<specimen>_<date>.txt`: MAPQ distributions in the padded window, soft-clipping near the junction and a test for reads counted twice (expected 0; also recorded in `run_provenance.txt`). To produce it for a `coverage_run/` folder from an older run:
  ```
  cd coverage_run
  python3 ../scripts/01_build_padded_reference.py --diagnose
  ```

## Outputs

When the driver runs `02_fold_and_plot.py`, the `<prefix>_` files below are in `coverage_results/`; run by hand, they are next to the script.


| File | Contents | What to look for |
| --- | --- | --- |
| `<prefix>_report.txt` | **Everything to report in the Methods and Results**: software versions, rotation and overlap checks, insert size, *K*, breadth, mean ± SD, median, minimum and maximum (with positions), both check results, and draft Methods and Results sentences built from them | Start here when writing |
| `<prefix>_coverage_stats.txt` | Genome length, *K*, mean, SD, minimum, maximum, breadth of coverage, end-window check result | Summary statistics |
| `<prefix>_end_window_check.txt` | Folded vs origin-check depth across the end windows (last and first *K* bp), in 200-bp bins | PASSED; if FAILED, the reported profile's end windows come from the origin-check mapping |
| `<prefix>_control_comparison.txt` | Corrected vs unpadded depth across the genome interior (positions *K*+1 to *N*−*K*) | Difference 0% or very close; same minimum at the same position(s) |
| `<prefix>_linear_coverage.png`, `<prefix>_circular_coverage.png` | Coverage figures (300 dpi) | No dip, spike or step at either end of the genome |
| `coverage_run/run_provenance.txt` | Date, software versions, settings, input checksums, check results (including both checks when the driver runs `02`) | Keep with the results |
| `coverage_run/diagnostic_report_<specimen>_<date>.txt` | MAPQ in the padded window, soft-clipping near the junction, double-count test | Reads present in both copies: 0; consult if a check fails |

The interior means in the comparison file leave out *K* bp at each end, so they differ from the whole-genome mean; report the whole-genome values from the statistics file.

## Validation

Padding changes only how reads near the origin are placed, so alignments in the genome interior are unaffected and the corrected and unpadded profiles should match there. In every specimen tested so far, interior mean depth agreed exactly or within 0.5%, with the same minimum at the same positions, across mean depths from about 40× to about 30,700×.

## End-window check

The interior comparison cannot see the first and last *K* bp, which is exactly where padding and folding act. In two *Hogna carolinensis* specimens rotated to begin at *cox1*, the folded profile showed excess depth over the first few hundred bases (at position 1, about 75% above the true value), although the interior comparison passed. Mapping the same reads to the unrotated assembly, where that region lies mid-sequence, gave a smooth profile. Simulated reads did not reproduce the effect, so it arises from features of particular read sets rather than from the folding arithmetic, and it cannot be ruled out in advance for any specimen.

Since version 4.2, `run_pipeline.sh` therefore also maps the reads (Phase 8), unpadded, to a copy of the assembly starting halfway round the circle, so the original start and end lie deep in its interior. `02_fold_and_plot.py` compares that profile with the folded profile across the end windows in 200-bp bins. A bin fails if it differs by 10% or more **and** by 3× or more. If no bin fails, the folded profile is reported as before. If any bin fails, the end windows of the reported profile are taken from the origin-check mapping, the statistics file says so, and the folded profile is kept as `<prefix>_folded_depth_before_end_window_correction.txt`. Either way, the reported profile stays in the assembly's own coordinates (for example *cox1* first).

## Provenance

Each run records how it was produced. `run_pipeline.sh` writes `run_provenance.txt` (software versions, settings, MD5 checksums of every input and script, check results, *K*, alignment counts). `padding_metadata.txt` records the assembly checksum, and the statistics and check files end with a provenance block giving the script version, date, Python library versions and input checksums. Mapping commands are also stored in each BAM header (`samtools view -H filtered.bam | grep "^@PG"`).

## Important notes

- **Do not add a MAPQ filter (`-q`).** Reads in the duplicated padding window align equally well to both copies and receive low mapping quality by design. Filtering them reintroduces the artifact. The proper-pair filter (`-f 2`) provides the quality control.
- **Use Bowtie2 in `--local` mode**, as the driver does. End-to-end mode is unforgiving at the padded junction.
- **Rotated assemblies:** the overlap check must be run on the assembler output before rotating, and the rotated file must be an exact cyclic shift. The driver requires an answer on whether the assembly was rotated and, if it was, runs both checks itself (`00_check_terminal_overlap.py`, with `--verify-rotation` for the second).

## Repository contents

| Path | Purpose |
| --- | --- |
| `scripts/run_pipeline.sh` | Driver: checks, mapping, padding, filtering, depth, diagnostic report, control and origin-check depth, provenance record; runs `02_fold_and_plot.py` if present |
| `scripts/00_check_terminal_overlap.py` | Checks for leftover circularization overlap and invalid characters; with `--verify-rotation`, confirms a rotated assembly is an exact cyclic shift of the original |
| `scripts/01_build_padded_reference.py` | Builds the padded reference and its metadata; with `--diagnose`, writes the padding diagnostic report |
| `scripts/02_fold_and_plot.py` | Folds depth, runs the end-window and interior checks, computes statistics, draws figures and writes the report; run by the driver (Phase 9) or by hand |
| `scripts/03_compare_to_control.py` | Optional standalone control comparison, for profiles folded with earlier versions |
| `docs/Coverage_pipeline_quick_start.docx` | Short quick start guide (four pages) |
| `docs/Coverage_calculation_and_visualization_protocol_v4.5.docx` | Full step-by-step protocol with embedded scripts |
| `docs/README_full_reference.txt` | Full reference: settings, outputs, troubleshooting, limitations |
| `docs/Supplement_1_command_line_protocol.docx` | Condensed command-line protocol for manuscript supplements |
| `docs/Coverage_Methods.docx` | Methods and Results template for manuscripts: a condensed Methods paragraph citing this repository (version 4.5), Results and figure-legend templates, and references |
| `CHANGELOG.md` | Version history |
| `CITATION.cff` | Citation information (powers GitHub's “Cite this repository” button) |

## Limitations

The pipeline expects a single circular contig and paired-end reads. The overlap check finds only exact suffix–prefix matches, so an overlap containing mismatches is not detected, and a genuine repeat at the origin can trigger it. Folding assumes one primary alignment per read, which holds with Bowtie2's default reporting (not with `-k`).

## How to cite this

If you use this pipeline, please cite it with the version you ran (recorded in `run_provenance.txt` and in the provenance block of each statistics file):

> Patrick LB. 2026. Mitogenome Coverage Pipeline, version 4.5. https://github.com/brpatric/mitogenome-coverage-pipeline

BibTeX:

```bibtex
@software{patrick_2026_mitogenome_coverage,
  author  = {Patrick, L. Brian},
  title   = {Mitogenome Coverage Pipeline},
  version = {4.5},
  year    = {2026},
  url     = {https://github.com/brpatric/mitogenome-coverage-pipeline}
}
```

GitHub's **Cite this repository** button (in the right-hand panel of the repository page) gives the same citation in APA and BibTeX formats, from `CITATION.cff`.

Please also cite the tools the pipeline runs: Bowtie2 (Langmead & Salzberg 2012), SAMtools (Danecek et al. 2021), NumPy (Harris et al. 2020), pandas (McKinney 2010) and Matplotlib (Hunter 2007); see References below. A ready-made Methods paragraph with these citations is in `docs/Coverage_Methods.docx`, and each run's `<prefix>_report.txt` gives draft Methods and Results text filled in with that run's values.

## References

1. Dierckxsens N, Mardulyn P, Smits G. 2017. NOVOPlasty: de novo assembly of organelle genomes from whole genome data. *Nucleic Acids Research* 45(4): e18. [doi:10.1093/nar/gkw955](https://doi.org/10.1093/nar/gkw955)
2. Langmead B, Salzberg SL. 2012. Fast gapped-read alignment with Bowtie 2. *Nature Methods* 9(4): 357–359. [doi:10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)
3. Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. 2021. Twelve years of SAMtools and BCFtools. *GigaScience* 10(2): giab008. [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)
4. Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, Cournapeau D, Wieser E, Taylor J, Berg S, Smith NJ, Kern R, Picus M, Hoyer S, van Kerkwijk MH, Brett M, Haldane A, Fernández del Río J, Wiebe M, Peterson P, Gérard-Marchant P, Sheppard K, Reddy T, Weckesser W, Abbasi H, Gohlke C, Oliphant TE. 2020. Array programming with NumPy. *Nature* 585(7825): 357–362. [doi:10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)
5. Hunter JD. 2007. Matplotlib: a 2D graphics environment. *Computing in Science & Engineering* 9(3): 90–95. [doi:10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55)
6. McKinney W. 2010. Data structures for statistical computing in Python. In: van der Walt S, Millman J, editors. *Proceedings of the 9th Python in Science Conference*. p. 56–61. [doi:10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a)
7. Jin JJ, Yu WB, Yang JB, Song Y, dePamphilis CW, Yi TS, Li DZ. 2020. GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes. *Genome Biology* 21(1): 241. [doi:10.1186/s13059-020-02154-5](https://doi.org/10.1186/s13059-020-02154-5)
