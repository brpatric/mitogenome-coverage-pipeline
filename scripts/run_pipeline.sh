#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh
#
# Circular-aware coverage mapping pipeline for a circularized
# mitogenome (developed with NOVOPlasty assemblies). Run this INSIDE WSL
# (Ubuntu). See README.txt for the full step-by-step explanation of what
# each phase does and why.
#
# Pipeline version 4.5. Phase 0 asks whether the assembly was rotated
# (or reads ROTATED below); if so, it runs the overlap check on the
# unrotated assembler output and verifies the rotation, and stops if either
# fails. Phase 6 writes the padding diagnostic report (01_build_padded_
# reference.py --diagnose). Phase 7 builds the unpadded control depth file
# (prelim_filtered_depth.txt) and Phase 8 maps the reads to a copy of the
# assembly starting halfway round the circle (origin_check_depth.txt), used
# by 02_fold_and_plot.py for the interior and end-window checks. If
# 02_fold_and_plot.py sits next to this script, Phase 9 runs it and writes
# the statistics, figures, check results and the report (everything to put
# in the Methods and Results) to RESULTS_DIR. Reads may
# be .fastq or .fastq.gz and are used in place (not copied); file settings
# accept names, relative paths or absolute paths. Writes run_provenance.txt
# (date, pipeline and tool versions, settings, input checksums and check
# results) to the output folder.
#
# Requires: bowtie2, samtools, python3 (all inside WSL)
# =============================================================================
set -euo pipefail

# --------------------------- USER-EDITABLE SETTINGS -------------------------
# File settings: a file name in the specimen folder (the folder you run
# this script from), a path relative to it, or an absolute path.
ASSEMBLY="assembly.fasta"      # your circularized mitogenome FASTA
R1="R1.fastq.gz"               # forward reads (.fastq or .fastq.gz)
R2="R2.fastq.gz"               # reverse reads (.fastq or .fastq.gz)
THREADS=8
OUTDIR="coverage_run"          # everything gets written here
MIN_PAD=300                    # floor on padding length K, in bp
MIN_OVERLAP=20                 # shortest terminal match treated as real
                               # redundancy (bp). Shorter matches occur by
                               # chance and are reported as a note only.
SD_MULTIPLIER=4                # K = mean_insert + SD_MULTIPLIER * sd_insert
ROTATED=""                     # "yes" if you rotated the assembly (for
                               # example to begin at cox1), "no" if it is the
                               # assembler output as is. Leave empty to be
                               # asked when the script starts; a run with no
                               # one to answer (e.g. a cluster job) stops
                               # instead of guessing.
RESULTS_DIR="coverage_results" # where 02_fold_and_plot.py's outputs go (a
                               # folder in the specimen folder). Used only if
                               # 02_fold_and_plot.py is next to this script.
PREFIX=""                      # output file prefix; empty = specimen folder
                               # name (e.g. Hogna_carolinensis_Mt19)
PLOT_TITLE="Per-base Read Coverage"  # title on both coverage figures
ORIGINAL_ASSEMBLY=""           # if ROTATED="yes": the unrotated assembler
                               # output (e.g. "original_assembly.fasta").
                               # Leave empty to be asked.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# -----------------------------------------------------------------------------
PIPELINE_VERSION="4.5"

# ------------------- rotated or not? (answer required) ----------------------
norm_yesno() { case "$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')" in
                   y|yes) echo yes ;; n|no) echo no ;; *) echo "" ;; esac; }
ROTATED_IN="$ROTATED"
ROTATED="$(norm_yesno "$ROTATED")"
if [ -n "$ROTATED_IN" ] && [ -z "$ROTATED" ]; then
    echo "!! ERROR: ROTATED is set to \"$ROTATED_IN\"; use \"yes\" or \"no\". Nothing was run."
    exit 1
fi
if [ -z "$ROTATED" ]; then
    if [ -t 0 ]; then
        while [ -z "$ROTATED" ]; do
            read -r -p "Was $ASSEMBLY rotated after assembly (e.g. to begin at cox1)? [yes/no]: " ans
            ROTATED="$(norm_yesno "$ans")"
            [ -z "$ROTATED" ] && echo "Please answer yes or no."
        done
    else
        echo "!! ERROR: ROTATED is not set and there is no terminal to ask."
        echo "!! Set ROTATED=\"yes\" or ROTATED=\"no\" at the top of run_pipeline.sh"
        echo "!! (and ORIGINAL_ASSEMBLY if yes). Nothing was run."
        exit 1
    fi
fi
if [ "$ROTATED" = "yes" ] && [ -z "$ORIGINAL_ASSEMBLY" ]; then
    if [ -t 0 ]; then
        while [ -z "$ORIGINAL_ASSEMBLY" ]; do
            read -r -p "Unrotated assembler output (file name or path): " ORIGINAL_ASSEMBLY
        done
    else
        echo "!! ERROR: ROTATED=\"yes\" but ORIGINAL_ASSEMBLY is not set, and there is"
        echo "!! no terminal to ask. Set ORIGINAL_ASSEMBLY at the top of"
        echo "!! run_pipeline.sh. Nothing was run."
        exit 1
    fi
fi
if [ "$ROTATED" = "no" ] && [ -n "$ORIGINAL_ASSEMBLY" ]; then
    echo "!! ERROR: ROTATED=\"no\" but ORIGINAL_ASSEMBLY is set. If the assembly was"
    echo "!! rotated, set ROTATED=\"yes\"; otherwise clear ORIGINAL_ASSEMBLY."
    echo "!! Nothing was run."
    exit 1
fi

# Resolve input paths against the specimen folder (where this is run from),
# before moving into the output folder. Absolute paths are kept as given.
abspath() { case "$1" in /*) printf '%s\n' "$1" ;; *) printf '%s\n' "$PWD/$1" ;; esac; }
ASSEMBLY_PATH="$(abspath "$ASSEMBLY")"
R1_PATH="$(abspath "$R1")"
R2_PATH="$(abspath "$R2")"
ORIGINAL_PATH=""
[ -n "$ORIGINAL_ASSEMBLY" ] && ORIGINAL_PATH="$(abspath "$ORIGINAL_ASSEMBLY")"

for f in "$ASSEMBLY_PATH" "$R1_PATH" "$R2_PATH" ${ORIGINAL_PATH:+"$ORIGINAL_PATH"}; do
    if [ ! -f "$f" ]; then
        echo "!! ERROR: input file not found: $f"
        echo "!! Check ASSEMBLY, R1, R2 and ORIGINAL_ASSEMBLY at the top of"
        echo "!! run_pipeline.sh. Nothing was run."
        exit 1
    fi
done

SPECIMEN_DIR="$PWD"
[ -z "$PREFIX" ] && PREFIX="$(basename "$SPECIMEN_DIR")"
RESULTS_PATH="$(abspath "$RESULTS_DIR")"
FOLD_SCRIPT="$SCRIPT_DIR/02_fold_and_plot.py"

mkdir -p "$OUTDIR"
cd "$OUTDIR"
# The assembly is small, so a working copy is kept with the results. The
# read files are used where they are (no copy), compressed or not.
cp "$ASSEMBLY_PATH" ./assembly.fasta

# ---------------------------- provenance record ------------------------------
PROV="run_provenance.txt"
{
echo "Circular-aware mitogenome coverage pipeline -- run provenance"
echo "Pipeline version   : $PIPELINE_VERSION"
echo "Run started        : $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "Working folder     : $(pwd)"
echo ""
echo "[Software]"
echo "bowtie2            : $(bowtie2 --version 2>/dev/null | head -1 | sed 's/.*version //')"
echo "samtools           : $(samtools --version 2>/dev/null | head -1 | sed 's/^samtools //')"
echo "python3            : $(python3 --version 2>&1 | sed 's/^Python //')"
echo ""
echo "[Settings]"
echo "ASSEMBLY           : $ASSEMBLY_PATH"
echo "ROTATED            : $ROTATED"
echo "ORIGINAL_ASSEMBLY  : ${ORIGINAL_PATH:-(not set)}"
echo "R1                 : $R1_PATH"
echo "R2                 : $R2_PATH"
echo "THREADS            : $THREADS"
echo "MIN_PAD            : $MIN_PAD"
echo "SD_MULTIPLIER      : $SD_MULTIPLIER"
echo "MIN_OVERLAP        : $MIN_OVERLAP"
if [ -f "$FOLD_SCRIPT" ]; then
echo "RESULTS_DIR        : $RESULTS_PATH"
echo "PREFIX             : $PREFIX"
echo "PLOT_TITLE         : $PLOT_TITLE"
fi
echo "bowtie2 mode       : --local (default scoring and reporting)"
echo "alignment filter   : samtools view -f 2 (properly paired), no MAPQ threshold"
echo ""
echo "[Input checksums (MD5)]"
for f in "$ASSEMBLY_PATH" "$R1_PATH" "$R2_PATH" ${ORIGINAL_PATH:+"$ORIGINAL_PATH"}; do
    printf "%s  %s\n" "$(md5sum "$f" | cut -d' ' -f1)" "$f"
done
echo ""
echo "[Scripts]"
for f in 00_check_terminal_overlap.py 01_build_padded_reference.py 02_fold_and_plot.py; do
    if [ -f "$SCRIPT_DIR/$f" ]; then
        printf "%-30s md5 %s\n" "$f" "$(md5sum "$SCRIPT_DIR/$f" | cut -d' ' -f1)"
    fi
done
printf "%-30s md5 %s\n" "run_pipeline.sh" "$(md5sum "$SCRIPT_DIR/run_pipeline.sh" | cut -d' ' -f1)"
echo ""
echo "[Results of checks and derived values]"
} > "$PROV"

echo "============================================================"
echo "PHASE 0: Checks on the assembly"
echo "============================================================"
OVERLAP_SCRIPT="$SCRIPT_DIR/00_check_terminal_overlap.py"
for f in "$OVERLAP_SCRIPT" "$SCRIPT_DIR/01_build_padded_reference.py"; do
    if [ ! -f "$f" ]; then
        echo "!! ERROR: required script not found:"
        echo "!!   $f"
        echo "!! It must sit in the same folder as run_pipeline.sh. Stopping here"
        echo "!! rather than silently skipping a check."
        exit 1
    fi
done

# overlap_check <fasta> <log> <label>: runs 00_check_terminal_overlap.py,
# records the result, and stops on redundancy or on a crash.
overlap_check() {
    set +e
    python3 "$OVERLAP_SCRIPT" "$1" --min-overlap "$MIN_OVERLAP" 2>&1 | tee "$2"
    local rc=${PIPESTATUS[0]}
    set -e
    # Exit 0 = no redundancy, 1 = redundancy found; anything else = crash.
    if [ "$rc" -ne 0 ] && [ "$rc" -ne 1 ]; then
        echo ""
        echo "!! ERROR: the overlap check on $1 did not run correctly (exit code $rc)."
        echo "!! See $2. Stopping rather than continuing with an unverified assembly."
        exit 1
    fi
    printf "%-21s: %s (%s)\n" "$3" "$(grep '^RESULT:' "$2" | sed 's/^RESULT: //')" "$2" >> "$PROV"
    if grep -q "RESULT: Exact terminal redundancy" "$2"; then
        echo ""
        echo "!! Terminal redundancy detected in $1 -- see $2."
        echo "!! Trim the reported overlap (from the assembler output, then redo"
        echo "!! any rotation) and re-run this pipeline. Stopping here."
        exit 1
    fi
}

if [ "$ROTATED" = "yes" ]; then
    cp "$ORIGINAL_PATH" ./original_assembly.fasta
    echo "--- Overlap check on the unrotated assembler output ---"
    overlap_check original_assembly.fasta overlap_check_original.log "Overlap (unrotated)"
    echo ""
    echo "--- Rotation check ---"
    set +e
    python3 "$OVERLAP_SCRIPT" --verify-rotation original_assembly.fasta assembly.fasta 2>&1 | tee rotation_check.log
    ROT_RC=${PIPESTATUS[0]}
    set -e
    if [ "$ROT_RC" -ne 0 ]; then
        echo "Rotation check       : FAILED (see rotation_check.log)" >> "$PROV"
        echo ""
        echo "!! The rotated assembly is NOT an exact cyclic shift of"
        echo "!! $ORIGINAL_ASSEMBLY. See rotation_check.log. Stopping."
        exit 1
    fi
    echo "Rotation check       : PASSED (rotation_check.log)" >> "$PROV"
    echo ""
    echo "--- Overlap check on the rotated assembly ---"
    overlap_check assembly.fasta overlap_check.log "Overlap check"
else
    echo "Rotation check       : not needed (assembly not rotated)" >> "$PROV"
    overlap_check assembly.fasta overlap_check.log "Overlap check"
fi

echo ""
echo "============================================================"
echo "PHASE 1: Preliminary mapping to estimate insert size"
echo "============================================================"
bowtie2-build assembly.fasta assembly_index
bowtie2 --local -p "$THREADS" -x assembly_index -1 "$R1_PATH" -2 "$R2_PATH" -S prelim.sam
samtools view -bS prelim.sam | samtools sort -o prelim.sorted.bam
samtools index prelim.sorted.bam
samtools stats prelim.sorted.bam > prelim.stats.txt

MEAN_INS=$(grep "^SN" prelim.stats.txt | grep "insert size average" | cut -f3)
SD_INS=$(grep "^SN" prelim.stats.txt | grep "insert size standard deviation" | cut -f3)

if [ -z "$MEAN_INS" ] || [ -z "$SD_INS" ]; then
    echo "Could not parse insert size from prelim.stats.txt -- falling back to MIN_PAD=$MIN_PAD"
    K=$MIN_PAD
else
    echo "Empirical insert size: mean=$MEAN_INS  sd=$SD_INS"
    K=$(python3 -c "
import math
mean = float('$MEAN_INS')
sd = float('$SD_INS')
k = math.ceil(mean + $SD_MULTIPLIER * sd)
k = max(k, $MIN_PAD)
print(k)
")
fi
echo "Selected padding length K = $K bp"
echo "$K" > padding_K.txt
{
echo "Insert size mean     : ${MEAN_INS:-(not parsed)}"
echo "Insert size SD       : ${SD_INS:-(not parsed)}"
echo "Padding length K     : $K bp"
} >> "$PROV"

echo ""
echo "============================================================"
echo "PHASE 2: Build padded reference"
echo "============================================================"
python3 "$SCRIPT_DIR/01_build_padded_reference.py" assembly.fasta "$K" padded_assembly.fasta padding_metadata.txt

echo ""
echo "============================================================"
echo "PHASE 3: Map reads to padded reference"
echo "============================================================"
bowtie2-build padded_assembly.fasta padded_assembly_index
bowtie2 --local -p "$THREADS" -x padded_assembly_index -1 "$R1_PATH" -2 "$R2_PATH" -S mapped.sam
samtools view -bS mapped.sam | samtools sort -o mapped.sorted.bam
samtools index mapped.sorted.bam

echo ""
echo "============================================================"
echo "PHASE 4: Filter for quality (proper pairs only)"
echo "============================================================"
samtools view -b -f 2 mapped.sorted.bam > filtered.bam
samtools index filtered.bam

echo ""
echo "============================================================"
echo "PHASE 5: Extract per-base depth from padded reference"
echo "============================================================"
samtools depth -a filtered.bam > padded_depth.txt

echo ""
echo "============================================================"
echo "PHASE 6: Padding diagnostic report"
echo "============================================================"
# 01_build_padded_reference.py --diagnose: MAPQ in the padded window,
# soft-clipping near the junction, and the double-count test, written to
# diagnostic_report_<specimen>_<date>.txt. A problem in this report-only
# step never stops the run.
set +e
python3 "$SCRIPT_DIR/01_build_padded_reference.py" --diagnose .
DIAG_RC=$?
set -e
DIAG_REPORT=$(ls -t diagnostic_report_*.txt 2>/dev/null | head -1)
if [ "$DIAG_RC" -eq 0 ] && [ -n "$DIAG_REPORT" ]; then
    DUPS=$(grep "reads present in BOTH copies" "$DIAG_REPORT" | sed 's/.*: *//')
    echo "Diagnostic report    : $DIAG_REPORT (reads in both copies of the window: ${DUPS:-?})" >> "$PROV"
else
    echo "!! WARNING: the diagnostic report did not complete (exit code $DIAG_RC)."
    echo "Diagnostic report    : did not complete (exit code $DIAG_RC)" >> "$PROV"
fi

echo ""
echo "============================================================"
echo "PHASE 7: Build unpadded control depth (for validation)"
echo "============================================================"
# The Phase 1 mapping to the UNPADDED assembly, filtered exactly as in
# Phase 4. 02_fold_and_plot.py compares it with the corrected profile
# across the genome interior (positions K+1 .. N-K).
samtools view -b -f 2 prelim.sorted.bam > prelim_filtered.bam
samtools index prelim_filtered.bam
samtools depth -a prelim_filtered.bam > prelim_filtered_depth.txt

N_EXPECTED=$(awk -F'\t' '$1=="original_length_N"{print $2}' padding_metadata.txt)
N_CONTROL=$(wc -l < prelim_filtered_depth.txt | tr -d '[:space:]')
if [ "$N_CONTROL" != "$N_EXPECTED" ]; then
    echo "!! WARNING: prelim_filtered_depth.txt has $N_CONTROL rows; expected N = $N_EXPECTED."
fi

echo ""
echo "============================================================"
echo "PHASE 8: Origin-check mapping (for the end-window check)"
echo "============================================================"
# The interior check (Phase 7) cannot see the ends of the sequence, which is
# exactly where padding acts. Here the reads are mapped, unpadded, to a copy
# of the assembly cut at the opposite side of the circle, so the original
# start/end region lies deep in its interior. 02_fold_and_plot.py compares
# this profile with the folded profile across the end windows.
S=$(python3 - <<'PY'
seq = "".join(l.strip() for l in open("assembly.fasta") if not l.startswith(">")).upper()
s = len(seq) // 2
with open("origin_check_assembly.fasta", "w") as fh:
    fh.write(f">origin_check_shift_{s}\n")
    rot = seq[s:] + seq[:s]
    for i in range(0, len(rot), 70):
        fh.write(rot[i:i + 70] + "\n")
print(s)
PY
)
echo "Origin-check assembly starts at original position $((S + 1))"
printf "origin_check_shift_S\t%s\n" "$S" >> padding_metadata.txt
bowtie2-build origin_check_assembly.fasta origin_check_index
bowtie2 --local -p "$THREADS" -x origin_check_index -1 "$R1_PATH" -2 "$R2_PATH" -S origin_check.sam
samtools view -bS origin_check.sam | samtools sort -o origin_check.sorted.bam
samtools view -b -f 2 origin_check.sorted.bam > origin_check_filtered.bam
samtools index origin_check_filtered.bam
samtools depth -a origin_check_filtered.bam > origin_check_depth.txt
N_ORIGIN=$(wc -l < origin_check_depth.txt | tr -d '[:space:]')
if [ "$N_ORIGIN" != "$N_EXPECTED" ]; then
    echo "!! WARNING: origin_check_depth.txt has $N_ORIGIN rows; expected N = $N_EXPECTED."
fi

{
echo "Genome length N      : $N_EXPECTED bp"
echo "Origin-check shift S : $S bp (origin-check assembly starts at position $((S + 1)))"
echo "Properly paired alignments (padded)   : $(samtools view -c filtered.bam)"
echo "Properly paired alignments (unpadded) : $(samtools view -c prelim_filtered.bam)"
echo "Properly paired alignments (origin check) : $(samtools view -c origin_check_filtered.bam)"
echo ""
echo "[Output checksums (MD5)]"
md5sum padded_depth.txt padding_metadata.txt prelim_filtered_depth.txt origin_check_depth.txt | awk '{printf "%-26s : %s\n", $2, $1}'
} >> "$PROV"

FOLD_RC=""
if [ -f "$FOLD_SCRIPT" ]; then
    echo ""
    echo "============================================================"
    echo "PHASE 9: Fold depth, run the checks, write statistics and figures"
    echo "============================================================"
    set +e
    MPLBACKEND=Agg python3 "$FOLD_SCRIPT" --input-dir . --output-dir "$RESULTS_PATH" \
        --prefix "$PREFIX" --title "$PLOT_TITLE"
    FOLD_RC=$?
    set -e
    STATS="$RESULTS_PATH/${PREFIX}_coverage_stats.txt"
    if [ "$FOLD_RC" -eq 0 ] && [ -f "$STATS" ]; then
        { echo ""; echo "[Results (Phase 9)]"; } >> "$PROV"
        echo "Results folder       : $RESULTS_PATH" >> "$PROV"
        echo "Report               : ${PREFIX}_report.txt" >> "$PROV"
        echo "End-window check     : $(grep '^End-window check:' "$STATS" | sed 's/^End-window check: //')" >> "$PROV"
        CMP="$RESULTS_PATH/${PREFIX}_control_comparison.txt"
        [ -f "$CMP" ] && echo "Interior check       : $(grep '^Difference' "$CMP" | sed 's/.*: //'); same minimum: $(grep '^Same minimum' "$CMP" | sed 's/.*: //')" >> "$PROV"
    else
        echo "!! WARNING: 02_fold_and_plot.py did not complete (exit code $FOLD_RC)."
        echo "!! The mapping outputs in $OUTDIR/ are complete; run 02_fold_and_plot.py"
        echo "!! by hand once the problem is fixed (often a missing Python package:"
        echo "!!   pip install numpy pandas matplotlib )."
        echo "Fold and plot (02)   : did not complete (exit code $FOLD_RC)" >> "$PROV"
    fi
fi

{
echo ""
echo "Run finished       : $(date '+%Y-%m-%d %H:%M:%S %Z')"
} >> "$PROV"

echo ""
echo "============================================================"
echo "DONE. Outputs written to: $OUTDIR/"
echo "  - padded_depth.txt            (raw per-base depth, length N+K)"
echo "  - padding_metadata.txt        (N, K -- needed by the fold/plot step)"
echo "  - prelim_filtered_depth.txt   (unpadded control depth, length N)"
echo "  - origin_check_depth.txt      (origin-check depth, length N)"
echo "  - overlap_check.log           (Phase 0 report; also"
echo "    overlap_check_original.log and rotation_check.log if rotated)"
echo "  - diagnostic_report_*.txt     (Phase 6 padding diagnostic report)"
echo "  - run_provenance.txt          (versions, settings, checksums, results)"
echo "  - prelim.stats.txt            (insert size stats used to pick K)"
echo ""
if [ "$FOLD_RC" = "0" ]; then
    echo "Statistics, figures and check results: $RESULTS_PATH/"
    echo "  - ${PREFIX}_report.txt              (everything to report, with draft text)"
    echo "  - ${PREFIX}_coverage_stats.txt      (summary statistics)"
    echo "  - ${PREFIX}_end_window_check.txt    ${PREFIX}_control_comparison.txt"
    echo "  - ${PREFIX}_linear_coverage.png     ${PREFIX}_circular_coverage.png"
else
    echo "Next: copy padded_depth.txt, padding_metadata.txt,"
    echo "prelim_filtered_depth.txt and origin_check_depth.txt next to"
    echo "02_fold_and_plot.py and run it (in WSL or in your Windows Python IDLE"
    echo "-- see README.txt). It folds the depth, runs the end-window and"
    echo "interior checks, and writes the statistics and figures. (Or put"
    echo "02_fold_and_plot.py next to run_pipeline.sh and the driver runs it.)"
fi
echo "============================================================"
