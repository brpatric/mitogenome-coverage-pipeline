================================================================================
CIRCULAR-AWARE MITOGENOME COVERAGE PIPELINE
================================================================================

Per-base coverage calculation and visualization for circularized mitochondrial
genomes assembled from paired-end Illumina data, with correction for the
coverage artifact introduced at the circular origin.

Version 4.2 (November 2026). See section 13 for the version history.


================================================================================
TABLE OF CONTENTS
================================================================================

   1. Overview
   2. The problem this pipeline solves
   3. Dependencies and versions
   4. Files in this package
   5. Input files you must supply
   6. How to run
        6.1  Part A - mapping and depth extraction (WSL)
        6.2  Part B - folding, plotting and validation
   7. Output files
   8. Validating the correction
   9. Troubleshooting
  10. Reporting in a manuscript
  11. Assumptions and limitations
  12. References
  13. Version history


================================================================================
1. OVERVIEW
================================================================================

This pipeline maps paired-end Illumina reads back to a circularized mitogenome
assembly and produces per-base depth statistics and publication-quality
coverage figures. It was developed for assemblies produced by NOVOPlasty [1],
but works with any correctly circularized single-contig assembly, including
GetOrganelle output.

The pipeline runs in two parts:

  Part A (WSL/Ubuntu)
      Indexing, mapping, filtering, and depth extraction, plus the unpadded
      control depth used for validation, driven by a single shell script.

  Part B (Windows Python/IDLE, or WSL)
      Folding the depth profile back onto true circular coordinates,
      computing statistics, plotting, and validating the correction against
      the unpadded control.

Typical runtime is a few minutes per specimen, dominated by read mapping.


================================================================================
2. THE PROBLEM THIS PIPELINE SOLVES
================================================================================

A circular mitogenome is a closed loop, but the FASTA file representing it is
linear: it has an artificial cut point where the last base meets the first.
Standard aligners treat that cut point as a genuine sequence end. Two things
follow:

  - Reads that span the biological origin cannot align contiguously.

  - Read PAIRS with one mate on either side of the origin appear, to the
    aligner, to be almost a full genome length apart, so they fail
    proper-pair filtering and are discarded entirely.

The result is an artificial coverage dip at both ends of the sequence, in a
region where nothing is biologically wrong. In a published figure this is
easily mistaken for a low-coverage or misassembled region.

THE CORRECTION
--------------
The reference is padded by appending the first K bp of the assembly to its own
end, where K is large enough to accommodate the longest expected fragment.
Reads spanning the origin then have a contiguous place to map. After depth
extraction, the depth in the appended region is added back onto the
corresponding positions at the start of the genome and the appended region
discarded, giving a continuous, artifact-free profile of the true genome.

WHY THE PADDING IS ONLY K BP AND NOT A FULL DUPLICATE
-----------------------------------------------------
Duplicating the entire genome would make every read in the genome ambiguous.
Padding by K confines that ambiguity to the small junction window, where it is
handled by the folding step.


================================================================================
3. DEPENDENCIES AND VERSIONS
================================================================================

INSIDE WSL (Ubuntu on Windows 11)
---------------------------------

  Tool            Version used   Purpose
  --------------  -------------  ------------------------------------------
  Bowtie2 [2]     v2.4.4         read alignment (--local mode)
  SAMtools [3]    v1.13          SAM/BAM handling, insert stats, depth
  Python          3.8+           helper scripts

  Install:

      sudo apt update
      sudo apt install bowtie2 samtools python3

  Record your actual versions for the manuscript:

      bowtie2 --version | head -1
      samtools --version | head -1

FOR THE VISUALIZATION STEP (WSL *or* Windows Python/IDLE)
---------------------------------------------------------

  Package         Purpose
  --------------  ------------------------------------------
  pandas [6]      reading the depth table
  NumPy [4]       depth arithmetic and statistics
  Matplotlib [5]  coverage figures

  Install:

      pip install pandas numpy matplotlib

No bioinformatics tools are needed for this step, which is why it can run
under Windows Python/IDLE rather than inside WSL.


================================================================================
4. FILES IN THIS PACKAGE
================================================================================

  run_pipeline.sh                  >>> RUN THIS (in WSL)
        Driver script. Runs Phases 0-7: checks, mapping, filtering, depth,
        the unpadded control depth and the origin-check depth.

  00_check_terminal_overlap.py     called automatically by run_pipeline.sh
        Verifies the assembly has no duplicated overlap at the circular
        join. Also validates sequence characters.

  01_build_padded_reference.py     called automatically by run_pipeline.sh
        Builds the padded reference and writes padding_metadata.txt.

  02_fold_and_plot.py              >>> RUN THIS (in IDLE or WSL)
        Folds depth to true coordinates, runs the end-window and interior
        checks, computes statistics and makes figures.

  verify_rotation.py               needed if you rotated the assembly
        Confirms the rotated assembly is an exact cyclic shift of the
        original assembler output. Run by the driver (Phase 0) when
        ORIGINAL_ASSEMBLY is set, or by hand (section 5).

  diagnose_padding_artifact.sh     optional (in WSL)
        Troubleshooting only: writes a diagnostic report on the padded
        window (MAPQ, soft-clipping, double-count test). See section 9.

  03_compare_to_control.py         optional (in IDLE or WSL)
        Standalone version of the control comparison now built into
        02_fold_and_plot.py. Kept for depth profiles folded with an earlier
        version of 02_fold_and_plot.py; it reports identical values.

  README.txt                       this file

IMPORTANT: 00_check_terminal_overlap.py and 01_build_padded_reference.py must
sit in the SAME FOLDER as run_pipeline.sh, as must verify_rotation.py if you
set ORIGINAL_ASSEMBLY. The driver stops with an explicit error if a required
script is missing, rather than silently skipping it.

Every script carries a version number (pipeline version 4.2), which it
prints or writes into its output; see section 7, PROVENANCE.


================================================================================
5. INPUT FILES YOU MUST SUPPLY
================================================================================

Place these in one working folder, with the scripts in a scripts/ subfolder:

    <specimen>/
        assembly.fasta        <- circularized mitogenome, single contig
        R1.fastq.gz           <- forward reads (.fastq or .fastq.gz)
        R2.fastq.gz           <- reverse reads (.fastq or .fastq.gz)
        scripts/
            run_pipeline.sh
            00_check_terminal_overlap.py
            01_build_padded_reference.py
            02_fold_and_plot.py
            03_compare_to_control.py
            verify_rotation.py
            diagnose_padding_artifact.sh

assembly.fasta REQUIREMENTS
---------------------------
  - A single circularized contig. Only the first FASTA record is used; a
    warning is printed if more are present.

  - Non-redundant: the last base should be immediately followed, circularly,
    by the first base, with no leftover overlap from the circularization
    step.

  - Only nucleotide or IUPAC ambiguity characters (ACGTUNRYSWKMBDHV).
    Anything else halts the run with the offending position reported.

IF YOU ROTATED THE ASSEMBLY (e.g., TO BEGIN AT cox1)
----------------------------------------------------
Do the checks in this order:

  1. Before rotating, check the assembler output for leftover overlap:

         python3 scripts/00_check_terminal_overlap.py original_assembly.fasta

  2. Rotate the assembly.

  3. Confirm the rotation is an exact cyclic shift of the original (same
     length, same composition, no bases duplicated, lost, or altered at the
     new start position):

         python3 scripts/verify_rotation.py \
             original_assembly.fasta assembly.fasta

     Alternatively, set ORIGINAL_ASSEMBLY in run_pipeline.sh (section
     6.1) and the driver runs this check in Phase 0, halts if it fails,
     and records the result in run_provenance.txt.

A cyclic shift of a non-redundant circle is itself non-redundant, so these
two checks together cover the rotated file. Phase 0 of the driver then
re-runs the overlap check on the rotated assembly.fasta; it should also
report no redundancy.


================================================================================
6. HOW TO RUN
================================================================================

--------------------------------------------------------------------------------
6.1  PART A - MAPPING AND DEPTH EXTRACTION (WSL)
--------------------------------------------------------------------------------

Open WSL, navigate to the specimen folder, and run:

    cd /mnt/c/path/to/<specimen>
    bash scripts/run_pipeline.sh

Windows drives are mounted under /mnt/c/... inside WSL.

Before the first run, open scripts/run_pipeline.sh and check the settings
block at the top:

  Setting            Default          Notes
  -----------------  ---------------  --------------------------------------
  ASSEMBLY           assembly.fasta   your mitogenome FASTA
  R1                 R1.fastq.gz      your forward reads
  R2                 R2.fastq.gz      your reverse reads
  THREADS            8                set to available CPU threads
  OUTDIR             coverage_run     all output lands here
  MIN_PAD            300              floor on padding length K
  SD_MULTIPLIER      4                K = mean insert + 4 x SD
  MIN_OVERLAP        20               shortest terminal match treated as
                                      real redundancy (bp)
  ORIGINAL_ASSEMBLY  (empty)          if you rotated the assembly, the
                                      unrotated assembler output; the
                                      driver then verifies the rotation

MIN_PAD and SD_MULTIPLIER rarely need changing.

FILE SETTINGS. Set ASSEMBLY, R1, R2 (and ORIGINAL_ASSEMBLY, if used) to
your own files. Each may be a file name in the specimen folder, a path
relative to it, or an absolute path, for example:

    R1="Sample1_R1.fastq.gz"
    R1="../shared_reads/Sample1_R1.fastq.gz"
    R1="/mnt/c/data/Sample1_R1.fastq.gz"

Reads may be plain (.fastq) or gzipped (.fastq.gz); they are used where
they are, not copied. The assembly is copied into coverage_run/ as
assembly.fasta. If any file is not found, the driver stops before running
anything and names the missing file.

WHAT THE DRIVER DOES
--------------------
  Phase 0   If ORIGINAL_ASSEMBLY is set, verifies that assembly.fasta is
            an exact rotation of it. Then checks the assembly for leftover
            circularization overlap and invalid characters. Halts if any
            check fails.

  Phase 1   Preliminary mapping to the unpadded assembly, purely to measure
            the empirical insert-size distribution.

  Phase 2   Calculates K from that distribution and builds the padded
            reference.

  Phase 3   Maps reads to the padded reference with Bowtie2 --local.

  Phase 4   Filters to properly paired alignments (-f 2).

  Phase 5   Extracts per-base depth across the padded reference.

  Phase 6   Filters the Phase 1 (unpadded) mapping exactly as in Phase 4
            and extracts its depth, giving the control used to validate the
            correction (prelim_filtered_depth.txt).

  Phase 7   Maps the reads, unpadded, to a copy of the assembly starting
            halfway round the circle, filters exactly as in Phase 4 and
            extracts its depth (origin_check_depth.txt), for the end-window
            check (section 8). The shift is recorded in padding_metadata.txt
            as origin_check_shift_S.

BOWTIE2 MUST BE RUN IN --local MODE
-----------------------------------
The default end-to-end mode requires the entire read to align, which is
unforgiving precisely at the padded junction; local mode soft-clips instead,
so reads there degrade gracefully.

NO MAPPING-QUALITY FILTER IS APPLIED, AND NONE SHOULD BE ADDED
--------------------------------------------------------------
The padded reference contains an intentional duplicate of the first K bp, so
reads landing in that window map equally well to either copy and correctly
receive low MAPQ. Filtering on MAPQ discards them and reintroduces exactly
the artifact the padding exists to remove. In one test specimen, a "-q 20"
threshold discarded 0.0% of interior alignments but 63-97% of alignments in
the padded window. The "-f 2" proper-pair filter provides the real quality
control here.

--------------------------------------------------------------------------------
6.2  PART B - FOLDING AND PLOTTING
--------------------------------------------------------------------------------

Copy these four files from coverage_run/ into the same folder as
02_fold_and_plot.py:

    padded_depth.txt
    padding_metadata.txt
    prelim_filtered_depth.txt   (needed for the interior check)
    origin_check_depth.txt      (needed for the end-window check)

Open the script and set:

    PLOT_TITLE = "Per-base Read Coverage: Genus species"

Then run it. No command-line arguments are needed:

    Windows Python/IDLE:  open the file and press F5 (Run > Run Module)
    WSL:                  python3 02_fold_and_plot.py

When prompted, type an output prefix for the specimen (for example
Genus_species_specimenID) and press Enter. All outputs are named with that
prefix and written next to the script.

If origin_check_depth.txt is present, the script runs the end-window check
(section 8) and writes <prefix>_end_window_check.txt. If
prelim_filtered_depth.txt is present, it also runs the interior check and
writes <prefix>_control_comparison.txt. If either file is absent, the script
prints a note, skips that check, and produces the other outputs as normal.



================================================================================
7. OUTPUT FILES
================================================================================

FROM PART A, IN coverage_run/
-----------------------------
  padded_depth.txt        Per-base depth across the padded reference
                          (N + K rows)
  padding_metadata.txt    N, K, and padded length - required by Part B
  padding_K.txt           The selected K, for the manuscript
  prelim.stats.txt        Insert-size statistics used to derive K
  run_provenance.txt      Run record: versions, settings, checksums,
                          check results (see PROVENANCE below)
  rotation_check.log      Rotation check report (only if
                          ORIGINAL_ASSEMBLY was set)
  overlap_check.log       Phase 0 report
  filtered.bam            Filtered alignments, retained for diagnostics
  prelim_filtered_depth.txt
                          Unpadded control depth (N rows) - used by Part B
                          for the interior check
  origin_check_depth.txt  Origin-check depth (N rows, in the shifted copy's
                          coordinates) - used by Part B for the end-window
                          check

FROM PART B
-----------
  <prefix>_folded_depth.txt       Reported per-base depth in the assembly's
                                  coordinates
  <prefix>_coverage_stats.txt     MEAN, SD, MIN, MAX, BREADTH - report these;
                                  also states the end-window check result
  <prefix>_linear_coverage.png    Linear coverage figure (300 dpi)
  <prefix>_circular_coverage.png  Circular/polar alternative
  <prefix>_control_comparison.txt Interior mean and minimum depth of the
                                  corrected and control profiles (only if
                                  prelim_filtered_depth.txt was present)
  <prefix>_end_window_check.txt   Folded vs origin-check depth across the end
                                  windows (only if origin_check_depth.txt
                                  was present)
  <prefix>_folded_depth_before_end_window_correction.txt
                                  The folded profile, kept for the record
                                  (only if the end-window check failed)

PROVENANCE
----------
Each run documents how it was produced:

  run_provenance.txt (Part A, coverage_run/)
      Date; pipeline version; Bowtie2, SAMtools and Python versions; all
      driver settings; MD5 checksums of the assembly, read files and
      scripts; rotation and overlap check results; insert-size mean and SD;
      K; N; the origin-check shift; numbers of properly paired alignments
      (padded, unpadded and origin check);
      and MD5 checksums of the depth and metadata files passed to Part B.

  padding_metadata.txt
      Also records the assembly file name, its MD5 checksum and the
      pipeline version.

  <prefix>_coverage_stats.txt, <prefix>_control_comparison.txt and
  <prefix>_end_window_check.txt (Part B)
      End with a "--- Provenance ---" block: producing script and pipeline
      version, date, Python, NumPy, pandas and Matplotlib versions, the
      assembly checksum (from padding_metadata.txt), and MD5 checksums of
      every input file. The values above the block are unchanged by it.

The checksums link each reported number back to the exact files it came
from: the "Input (MD5)" lines in the statistics file match the "Output
checksums" in run_provenance.txt for the same run. Mapping commands and tool
versions are also stored in each BAM header:

    samtools view -H filtered.bam | grep "^@PG"

DO NOT COMPUTE STATISTICS FROM padded_depth.txt. It still contains the padded
region and will give skewed values. All reported statistics must come from
<prefix>_coverage_stats.txt, which is derived from the reported profile.


================================================================================
8. VALIDATING THE CORRECTION
================================================================================

The padding-and-folding approach can be checked against a control mapping to
the unpadded assembly using identical reads and filtering. Across the genome
interior - the region padding does not touch - the two approaches should
agree closely.

The comparison is restricted to positions K+1 to N-K, because the unpadded
control is depressed at both ends by the origin artifact. It runs
automatically: Phase 6 of the driver builds the control, and
02_fold_and_plot.py compares it with the corrected profile when
prelim_filtered_depth.txt sits next to the script. In
<prefix>_control_comparison.txt:

  Difference (corrected - control)   percent difference in interior mean
                                     depth; expected to be 0% or very close
  Same minimum at same position(s)   expected YES

Padding changes only how reads near the origin are placed, so alignments in
the interior are normally identical and the two profiles match exactly. In
every specimen tested so far, agreement was exact or within 0.5%, with the
same minimum at the same positions, across mean depths from about 40x to
about 30,700x. A larger difference indicates that the two mappings were not
run on the same reads or with the same filter.

Note that the interior means in this file are NOT the whole-genome means in
<prefix>_coverage_stats.txt. Report whole-genome statistics in a manuscript.

THE END-WINDOW CHECK (version 4.2+)
-----------------------------------
The interior comparison cannot see the first and last K bp, which is exactly
where padding and folding act. In two Hogna carolinensis specimens assembled
with NOVOPlasty and rotated to begin at cox1, the folded profile showed
excess depth over the first few hundred bases (up to 50% above the true
value at position 1), although the interior comparison passed. Mapping the
same reads to the unrotated assembly, where that region lies mid-sequence,
gave a smooth profile. Simulated reads did not reproduce the effect, so it
arises from features of particular read sets rather than from the folding
arithmetic, and it cannot be ruled out in advance.

Phase 7 of the driver therefore maps the reads, unpadded, to a copy of the
assembly starting halfway round the circle, where the original start and
end lie deep in the interior. 02_fold_and_plot.py compares that profile with
the folded profile across the end windows (positions N-W+1 to N and 1 to W,
where W = K, reduced if needed to stay 500 bp from the copy's own ends), in
200-bp bins. In <prefix>_end_window_check.txt:

  Bins failing   a bin fails if it differs by >= 10% AND >= 3x
  RESULT         PASSED: the folded profile is reported
                 FAILED: the end windows of the reported profile are taken
                 from the origin-check mapping; the folded profile is kept as
                 <prefix>_folded_depth_before_end_window_correction.txt

The file also gives the mean depth of both profiles across the windows, the
largest one-base difference, and the step in depth between position N and
position 1 in each profile. The thresholds are set near the top of
02_fold_and_plot.py (END_MAX_PCT, END_MIN_ABS, END_BIN_BP). Bins inside
long tandem repeats can differ slightly between the two mappings, because
reads there are placed at random among the copies; a failure confined to
such a region is harmless, since the reported profile then simply uses the
origin-check values there. The reported profile is always in the assembly's
own coordinates.


================================================================================
9. TROUBLESHOOTING
================================================================================

ERROR: "padded_depth has X rows but metadata says N+K = Y"
----------------------------------------------------------
The depth file and metadata describe different references. Usually one is
stale from an earlier run. Check that padding_metadata.txt,
padded_assembly.fasta, and the BAM header agree:

    grep -v "^>" padded_assembly.fasta | tr -d '\n' | wc -c
    samtools view -H filtered.bam | grep "^@SQ"
    cat padding_metadata.txt

All three should give the same padded length. Re-run from the first step
where they diverge. This error also appears if the assembly contained
characters that Python counted but Bowtie2 rejected - see below.

RUN HALTS REPORTING NON-NUCLEOTIDE CHARACTERS
---------------------------------------------
The assembly contains something outside ACGTUNRYSWKMBDHV. The reported
positions tell you where. These characters must be resolved before mapping:
Bowtie2 rejects them, silently shortening the indexed reference. Note that
some assemblers use non-nucleotide characters deliberately as annotations
rather than as sequence, so consult your assembler's documentation before
deciding how to handle them - deleting a flagged position shifts all
downstream coordinates.

overlap_check.log IS EMPTY, OR PHASE 0 REPORTS A MISSING SCRIPT
---------------------------------------------------------------
00_check_terminal_overlap.py is not in the same folder as run_pipeline.sh.
Earlier versions of the driver skipped the check silently in this situation;
the current version halts. Restore the file and re-run.

"bowtie2-build: command not found"
----------------------------------
The pipeline is being run under Windows rather than inside WSL, or Bowtie2 is
not installed. Part A must run inside WSL.

END-WINDOW CHECK FAILED
-----------------------
No action is needed to obtain correct values: the reported profile already
takes the end windows from the origin-check mapping, and the statistics file
says so. Look at the failing bins in <prefix>_end_window_check.txt. Failures
near position 1 or near position K point to the padded junction; failures
confined to a long tandem repeat reflect random placement of reads among the
repeat copies. To investigate the junction, run diagnose_padding_artifact.sh
(below) or view filtered.bam in a genome browser across positions N-200 to
N+200 of the padded reference.

COVERAGE LOOKS DEPRESSED NEAR POSITION 1
----------------------------------------
Check <prefix>_end_window_check.txt first. Confirm no "-q" MAPQ threshold has been added to the filtering step (see
section 6.1). If the filter is clean, check <prefix>_control_comparison.txt,
then run diagnose_padding_artifact.sh from inside the specimen's
coverage_run/ folder and inspect its report (MAPQ in the padded window,
soft-clipping near the junction, and the double-count test):

    cd /mnt/c/path/to/<specimen>/coverage_run
    cp ../scripts/diagnose_padding_artifact.sh .
    bash diagnose_padding_artifact.sh

02_fold_and_plot.py REPORTS MISSING INPUT FILES
-----------------------------------------------
It looks for padded_depth.txt and padding_metadata.txt in the folder
containing the script itself, not the current working directory. Copy both
files next to the script.


================================================================================
10. REPORTING IN A MANUSCRIPT
================================================================================

A Methods/Results template (Coverage_Methods.docx) accompanies this
package. Keep the procedure in the Methods and the numbers in the Results.

METHODS (procedure only)
  - Overlap check before rotation, and rotation verification, if the
    assembly was rotated
  - Aligner and version (Bowtie2 v2.4.4, --local mode) [2] and SAMtools
    version [3]
  - How K was set (mean insert size + 4 SD, minimum 300 bp), padding, and
    folding of depth back to true coordinates before statistics
  - Filtering: properly paired reads only, NO MAPQ THRESHOLD, with the
    reason (section 6.1). An omitted quality filter looks like an oversight
    unless the rationale is given.
  - Analysis libraries: NumPy [4], pandas [6], Matplotlib [5]. Software
    versions are in run_provenance.txt (Bowtie2, SAMtools) and in the
    provenance block of <prefix>_coverage_stats.txt (Python libraries).
  - How the correction was validated: interior check and end-window check
    (section 8)
  - A pointer to the supplement containing the commands and scripts

RESULTS (values)
  - K for each specimen                 <prefix>_coverage_stats.txt
  - Mean +/- SD, range, breadth         <prefix>_coverage_stats.txt
  - Validation outcome (interior)       <prefix>_control_comparison.txt
  - End-window check result            <prefix>_end_window_check.txt
                                       (and, if it failed, that the end
                                       windows came from the origin-check
                                       mapping)

Padding and folding is less familiar than plain linear mapping, so reporting
the validation outcome pre-empts questions about whether the correction
introduced its own bias. If a few high-depth peaks make the SD larger than
the mean, consider also reporting the median depth (from
<prefix>_folded_depth.txt).


================================================================================
11. ASSUMPTIONS AND LIMITATIONS
================================================================================

  - SINGLE CIRCULAR CONTIG. Only the first FASTA record is processed. Linear
    or fragmented assemblies are out of scope; for a linear genome, no
    padding is needed and standard mapping suffices.

  - NON-REDUNDANT CIRCLE. The pipeline assumes the assembly does not already
    duplicate sequence at the join. Phase 0 checks this by exact
    suffix/prefix match, which will not detect an inexact
    (mismatch-containing) overlap. Only matches of at least MIN_OVERLAP bp
    (default 20) are treated as redundancy; shorter matches occur by chance
    - the last base of any sequence matches the first about one time in
    four - and are reported as a note only. Adjust MIN_OVERLAP in
    run_pipeline.sh if your assembler leaves shorter overlaps.

  - TERMINAL REPEATS. Genomes with genuine long repeats at the assembly
    origin may trigger the Phase 0 warning spuriously; inspect before
    trimming.

  - INSERT-SIZE-DERIVED K. K is set from the empirical insert-size
    distribution. A strongly skewed or bimodal library may need a larger
    SD_MULTIPLIER; check prelim.stats.txt if coverage near the origin looks
    irregular.

  - ONE PRIMARY ALIGNMENT PER READ. Folding assumes each read contributes
    depth once. This holds under the default Bowtie2 settings used here.
    Enabling multi-mapping reporting (-k) would break it.

  - PAIRED-END DATA. The -f 2 proper-pair filter assumes paired reads.
    Single-end data requires a different filtering strategy.


================================================================================
12. REFERENCES
================================================================================

[1] Dierckxsens N, Mardulyn P, Smits G. 2017. NOVOPlasty: de novo assembly of
    organelle genomes from whole genome data. Nucleic Acids Research
    45(4):e18. doi:10.1093/nar/gkw955

[2] Langmead B, Salzberg SL. 2012. Fast gapped-read alignment with Bowtie 2.
    Nature Methods 9(4):357-359. doi:10.1038/nmeth.1923

[3] Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO,
    Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. 2021. Twelve years of
    SAMtools and BCFtools. GigaScience 10(2):giab008.
    doi:10.1093/gigascience/giab008

[4] Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P,
    Cournapeau D, Wieser E, Taylor J, Berg S, Smith NJ, Kern R, Picus M,
    Hoyer S, van Kerkwijk MH, Brett M, Haldane A, Fernandez del Rio J,
    Wiebe M, Peterson P, Gerard-Marchant P, Sheppard K, Reddy T,
    Weckesser W, Abbasi H, Gohlke C, Oliphant TE. 2020. Array programming
    with NumPy. Nature 585(7825):357-362. doi:10.1038/s41586-020-2649-2

[5] Hunter JD. 2007. Matplotlib: a 2D graphics environment. Computing in
    Science & Engineering 9(3):90-95. doi:10.1109/MCSE.2007.55

[6] McKinney W. 2010. Data structures for statistical computing in Python.
    Proceedings of the 9th Python in Science Conference, p. 56-61.
    doi:10.25080/Majora-92bf1922-00a


================================================================================
13. VERSION HISTORY
================================================================================

  Version 4.2 (November 2026)
    - End-window check: Phase 7 of run_pipeline.sh maps the reads to a copy
      of the assembly starting halfway round the circle
      (origin_check_depth.txt); 02_fold_and_plot.py compares it with the
      folded profile across the first and last K bp and, if they disagree,
      takes the end windows from it (section 8).
    - The statistics file states the end-window check result.
    - Profiles that pass the check are unchanged from version 4.1.1.

  Version 4.1.1 (November 2026)
    - Reads may be .fastq or .fastq.gz and are used in place (not copied
      into coverage_run/).
    - File settings accept names, relative paths or absolute paths; a
      missing file stops the run with a clear message.
    - run_provenance.txt records full input paths.
    - Computed values are unchanged from version 4.1.

  Version 4.1 (November 2026)
    - Provenance: run_pipeline.sh writes run_provenance.txt; every script
      carries a version number; padding_metadata.txt records the assembly
      checksum; statistics and comparison files end with a provenance
      block (versions, date, input checksums).
    - Optional ORIGINAL_ASSEMBLY setting: the driver verifies a rotated
      assembly in Phase 0 and records the result.
    - Computed values (K, folded depth, statistics, comparison) are
      unchanged from version 4.

  Version 4 (November 2026)
    - run_pipeline.sh Phase 6 builds the unpadded control depth
      (prelim_filtered_depth.txt) automatically.
    - 02_fold_and_plot.py runs the control comparison when that file is
      present; output identical to 03_compare_to_control.py, which is kept
      as an optional standalone script.
    - diagnose_padding_artifact.sh is now a troubleshooting tool only.
    - Folded depth, statistics and figures are unchanged from version 3.

  Version 3 (November 2026)
    - Added 03_compare_to_control.py and Part C (section 6.3) to validate
      the correction against an unpadded control.
    - Overlap check now run on the assembler output before rotation,
      followed by verify_rotation.py.
    - Reporting guidance split into Methods (procedure) and Results
      (values); Methods/Results template and Supplement 1 added.
    - Bowtie2 version recorded as v2.4.4; full author list for [4].
    - Documentation made generic (no specimen-specific examples).

  Version 2
    - Folding of padded depth, terminal-overlap check, character
      validation, MAPQ guidance, verify_rotation.py and
      diagnose_padding_artifact.sh.


================================================================================
END OF README
================================================================================
