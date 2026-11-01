#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh
#
# Circular-aware coverage mapping pipeline for a circularized
# mitogenome (developed with NOVOPlasty assemblies). Run this INSIDE WSL
# (Ubuntu). See README.txt for the full step-by-step explanation of what
# each phase does and why.
#
# Pipeline version 4.1. Writes run_provenance.txt (date, pipeline and tool
# versions, settings, input checksums and check results) to the output
# folder. Phase 6 builds the unpadded control depth file
# (prelim_filtered_depth.txt) used for validation.
#
# Requires: bowtie2, samtools, python3 (all inside WSL)
# =============================================================================
set -euo pipefail

# --------------------------- USER-EDITABLE SETTINGS -------------------------
ASSEMBLY="assembly.fasta"      # your circularized mitogenome FASTA
R1="R1.fastq"                  # forward reads
R2="R2.fastq"                  # reverse reads
THREADS=8
OUTDIR="coverage_run"          # everything gets written here
MIN_PAD=300                    # floor on padding length K, in bp
MIN_OVERLAP=20                 # shortest terminal match treated as real
                               # redundancy (bp). Shorter matches occur by
                               # chance and are reported as a note only.
SD_MULTIPLIER=4                # K = mean_insert + SD_MULTIPLIER * sd_insert
ORIGINAL_ASSEMBLY=""           # OPTIONAL: if you rotated the assembly, the
                               # unrotated assembler output (e.g.
                               # "original_assembly.fasta"). When set, the
                               # driver verifies the rotation (Phase 0) and
                               # records the result. Leave empty otherwise.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# -----------------------------------------------------------------------------
PIPELINE_VERSION="4.1"

mkdir -p "$OUTDIR"
cd "$OUTDIR"
cp "../$ASSEMBLY" ./assembly.fasta
cp "../$R1" ./R1.fastq
cp "../$R2" ./R2.fastq

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
echo "ASSEMBLY           : $ASSEMBLY"
echo "ORIGINAL_ASSEMBLY  : ${ORIGINAL_ASSEMBLY:-(not set)}"
echo "R1                 : $R1"
echo "R2                 : $R2"
echo "THREADS            : $THREADS"
echo "MIN_PAD            : $MIN_PAD"
echo "SD_MULTIPLIER      : $SD_MULTIPLIER"
echo "MIN_OVERLAP        : $MIN_OVERLAP"
echo "bowtie2 mode       : --local (default scoring and reporting)"
echo "alignment filter   : samtools view -f 2 (properly paired), no MAPQ threshold"
echo ""
echo "[Input checksums (MD5)]"
md5sum assembly.fasta R1.fastq R2.fastq | awk '{printf "%-18s : %s\n", $2, $1}'
if [ -n "$ORIGINAL_ASSEMBLY" ]; then
    printf "%-18s : %s\n" "$ORIGINAL_ASSEMBLY" "$(md5sum "../$ORIGINAL_ASSEMBLY" | cut -d' ' -f1)"
fi
echo ""
echo "[Scripts]"
for f in 00_check_terminal_overlap.py 01_build_padded_reference.py verify_rotation.py; do
    if [ -f "$SCRIPT_DIR/$f" ]; then
        printf "%-30s md5 %s\n" "$f" "$(md5sum "$SCRIPT_DIR/$f" | cut -d' ' -f1)"
    fi
done
printf "%-30s md5 %s\n" "run_pipeline.sh" "$(md5sum "$SCRIPT_DIR/run_pipeline.sh" | cut -d' ' -f1)"
echo ""
echo "[Results of checks and derived values]"
} > "$PROV"

echo "============================================================"
echo "PHASE 0: Check for leftover terminal redundancy"
echo "============================================================"
OVERLAP_SCRIPT="$SCRIPT_DIR/00_check_terminal_overlap.py"
if [ ! -f "$OVERLAP_SCRIPT" ]; then
    echo "!! ERROR: required script not found:"
    echo "!!   $OVERLAP_SCRIPT"
    echo "!! 00_check_terminal_overlap.py must sit in the same folder as"
    echo "!! run_pipeline.sh. Stopping here rather than silently skipping"
    echo "!! the terminal-redundancy check."
    exit 1
fi

if [ -n "$ORIGINAL_ASSEMBLY" ]; then
    ROT_SCRIPT="$SCRIPT_DIR/verify_rotation.py"
    if [ ! -f "$ROT_SCRIPT" ]; then
        echo "!! ERROR: ORIGINAL_ASSEMBLY is set but verify_rotation.py was not"
        echo "!! found in $SCRIPT_DIR. Stopping."
        exit 1
    fi
    cp "../$ORIGINAL_ASSEMBLY" ./original_assembly.fasta
    set +e
    python3 "$ROT_SCRIPT" original_assembly.fasta assembly.fasta 2>&1 | tee rotation_check.log
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
else
    echo "Rotation check       : not run by the driver (ORIGINAL_ASSEMBLY not set)" >> "$PROV"
fi

set +e
python3 "$OVERLAP_SCRIPT" assembly.fasta --min-overlap "$MIN_OVERLAP" 2>&1 | tee overlap_check.log
OVERLAP_RC=${PIPESTATUS[0]}
set -e

# The script exits 0 (no redundancy) or 1 (redundancy found). Anything
# else means it crashed -- do not treat that as "all clear".
if [ "$OVERLAP_RC" -ne 0 ] && [ "$OVERLAP_RC" -ne 1 ]; then
    echo ""
    echo "!! ERROR: the overlap check did not run correctly (exit code $OVERLAP_RC)."
    echo "!! See overlap_check.log above. Stopping rather than continuing"
    echo "!! with an unverified assembly."
    exit 1
fi

echo "Overlap check        : $(grep '^RESULT:' overlap_check.log | sed 's/^RESULT: //') (overlap_check.log)" >> "$PROV"

if grep -q "RESULT: Exact terminal redundancy" overlap_check.log; then
    echo ""
    echo "!! Terminal redundancy detected -- see overlap_check.log for details."
    echo "!! Trim the reported overlap from assembly.fasta and re-run this"
    echo "!! pipeline before trusting the results. Stopping here."
    exit 1
fi

echo ""
echo "============================================================"
echo "PHASE 1: Preliminary mapping to estimate insert size"
echo "============================================================"
bowtie2-build assembly.fasta assembly_index
bowtie2 --local -p "$THREADS" -x assembly_index -1 R1.fastq -2 R2.fastq -S prelim.sam
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
bowtie2 --local -p "$THREADS" -x padded_assembly_index -1 R1.fastq -2 R2.fastq -S mapped.sam
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
echo "PHASE 6: Build unpadded control depth (for validation)"
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

{
echo "Genome length N      : $N_EXPECTED bp"
echo "Properly paired alignments (padded)   : $(samtools view -c filtered.bam)"
echo "Properly paired alignments (unpadded) : $(samtools view -c prelim_filtered.bam)"
echo ""
echo "[Output checksums (MD5)]"
md5sum padded_depth.txt padding_metadata.txt prelim_filtered_depth.txt | awk '{printf "%-26s : %s\n", $2, $1}'
echo ""
echo "Run finished       : $(date '+%Y-%m-%d %H:%M:%S %Z')"
} >> "$PROV"

echo ""
echo "============================================================"
echo "DONE. Outputs written to: $OUTDIR/"
echo "  - padded_depth.txt            (raw per-base depth, length N+K)"
echo "  - padding_metadata.txt        (N, K -- needed by the fold/plot step)"
echo "  - prelim_filtered_depth.txt   (unpadded control depth, length N)"
echo "  - overlap_check.log           (Phase 0 report)"
echo "  - run_provenance.txt          (versions, settings, checksums, results)"
echo "  - prelim.stats.txt            (insert size stats used to pick K)"
echo ""
echo "Next: copy padded_depth.txt, padding_metadata.txt and"
echo "prelim_filtered_depth.txt next to 02_fold_and_plot.py and run it"
echo "(in WSL or in your Windows Python IDLE -- see README.txt). It folds"
echo "the depth, writes the statistics and figures, and runs the control"
echo "comparison."
echo "============================================================"
