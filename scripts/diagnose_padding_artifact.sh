#!/usr/bin/env bash
# =============================================================================
# diagnose_padding_artifact.sh
#
# OPTIONAL troubleshooting tool. Use it if the finished coverage profile
# shows an unexplained dip or spike within the first K bases. It reports
# MAPQ distributions, soft-clipping and a double-count test for the padded
# window.
#
# Since pipeline version 4, run_pipeline.sh builds the unpadded control
# depth file (prelim_filtered_depth.txt) itself, so this script is no longer
# needed for routine validation. It rebuilds that file anyway (same
# commands, same result).
#
# HOW TO RUN
# ----------
#   1. Copy this file into a specimen's  coverage_run/  folder
#      (the folder containing padded_depth.txt, filtered.bam, etc.)
#   2. In WSL/Ubuntu:
#          cd /path/to/<specimen>/coverage_run
#          bash diagnose_padding_artifact.sh
#
# It does NOT modify or delete any existing pipeline output. It only reads
# existing files and writes two new ones (listed at the end).
#
# Requires: samtools (already needed by the main pipeline)
# =============================================================================
set -uo pipefail
PIPELINE_VERSION="4.1.1"

SPECIMEN="$(basename "$(dirname "$(pwd)")")"
STAMP="$(date +%Y%m%d)"
REPORT="diagnostic_report_${SPECIMEN}_${STAMP}.txt"

# ----------------------------------------------------------------------------
# Pre-flight: make sure we are in the right folder and have what we need
# ----------------------------------------------------------------------------
missing=0
for f in padded_depth.txt padding_metadata.txt filtered.bam prelim.sorted.bam; do
    if [ ! -f "$f" ]; then
        echo "MISSING: $f"
        missing=1
    fi
done
if [ "$missing" -ne 0 ]; then
    echo ""
    echo "This script must be run from inside a specimen's coverage_run/ folder."
    echo "Current folder: $(pwd)"
    echo "Nothing was written. Move to the correct folder and re-run."
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools not found on PATH. Run this inside WSL, not Windows."
    exit 1
fi

# K: prefer padding_K.txt, else parse it out of padding_metadata.txt
if [ -f padding_K.txt ]; then
    K="$(tr -d '[:space:]' < padding_K.txt)"
else
    K="$(awk -F'\t' '$1=="padding_length_K"{print $2}' padding_metadata.txt | tr -d '[:space:]')"
fi
TOTAL="$(wc -l < padded_depth.txt | tr -d '[:space:]')"
N=$((TOTAL - K))

# ----------------------------------------------------------------------------
# Everything from here is written to screen AND to the report file
# ----------------------------------------------------------------------------
{
echo "============================================================"
echo "PADDING ARTIFACT DIAGNOSTIC (pipeline version $PIPELINE_VERSION)"
echo "Specimen folder : $SPECIMEN"
echo "Path            : $(pwd)"
echo "Date            : $(date)"
echo "samtools        : $(samtools --version 2>/dev/null | head -1)"
echo "============================================================"
echo ""

echo "--- geometry ---"
echo "padded_depth.txt rows (N+K) : $TOTAL"
echo "K (padding length)          : $K"
echo "N (true genome length)      : $N"
echo ""
echo "padding_metadata.txt contents:"
sed 's/^/    /' padding_metadata.txt
echo ""

echo "--- BAM header reference length (should equal N+K = $TOTAL) ---"
samtools view -H filtered.bam | grep '^@SQ' | sed 's/^/    /'
echo ""

echo "--- alignment record counts in filtered.bam ---"
echo "total records          : $(samtools view -c filtered.bam)"
echo "secondary (flag 256)   : $(samtools view -c -f 256 filtered.bam)"
echo "supplementary (2048)   : $(samtools view -c -f 2048 filtered.bam)"
echo "duplicates (1024)      : $(samtools view -c -f 1024 filtered.bam)"
echo ""

echo "--- DOUBLE-COUNT TEST ---"
echo "Individual READS (QNAME + mate number) with alignments in BOTH copies"
echo "of the duplicated window (copy 1 = pos 1..$K, copy 2 = pos > $N)."
echo "NOTE: keying on QNAME alone is WRONG here -- both mates of a pair share"
echo "a QNAME, so an ordinary proper pair with both mates inside the window"
echo "would be counted as a false positive. Mate number is included below."
echo "If secondary/supplementary counts above are 0, this MUST be 0: with one"
echo "alignment record per read, a read cannot contribute depth twice."
DUPREADS=$(samtools view filtered.bam \
    | awk -v k="$K" -v n="$N" '
        { m1=int($2/64)%2; m2=int($2/128)%2;
          mate = m1 ? "/1" : (m2 ? "/2" : "/0");
          if ($4<=k)      print $1 mate "\tA";
          else if ($4>n)  print $1 mate "\tB" }' \
    | sort -u | cut -f1 | sort | uniq -d | wc -l)
echo "reads present in BOTH copies : $DUPREADS"
echo ""

echo "--- MAPQ distribution: window copy 1 (pos 1..$K) ---"
echo "   count  MAPQ"
samtools view filtered.bam | awk -v k="$K" '$4<=k {print $5}' \
    | sort -n | uniq -c | sed 's/^/    /'
echo ""

echo "--- MAPQ distribution: window copy 2 (pos > $N) ---"
echo "   count  MAPQ"
samtools view filtered.bam | awk -v n="$N" '$4>n {print $5}' \
    | sort -n | uniq -c | sed 's/^/    /'
echo ""

echo "--- MAPQ distribution: interior (pos $K..$N, first 200k records) ---"
echo "   count  MAPQ"
samtools view filtered.bam | awk -v k="$K" -v n="$N" '$4>k && $4<=n {print $5}' \
    | head -200000 | sort -n | uniq -c | sed 's/^/    /'
echo ""

echo "--- CIGAR soft-clipping near the junction (pos > $((N-200))) ---"
echo "   count  CIGAR (top 15)"
samtools view filtered.bam | awk -v x="$((N-200))" '$4>x {print $6}' \
    | sort | uniq -c | sort -rn | head -15 | sed 's/^/    /'
echo ""

echo "============================================================"
echo "END OF DIAGNOSTIC"
echo "============================================================"
} 2>&1 | tee "$REPORT"

# ----------------------------------------------------------------------------
# Unpadded control depth file -- the decisive comparison
# ----------------------------------------------------------------------------
echo ""
echo "Building unpadded control depth file (this may take a minute)..."
samtools view -b -f 2 prelim.sorted.bam > prelim_filtered.bam
samtools index prelim_filtered.bam
samtools depth -a prelim_filtered.bam > prelim_filtered_depth.txt
echo "  rows in prelim_filtered_depth.txt : $(wc -l < prelim_filtered_depth.txt)"
echo "  (should be N = $N, i.e. the UNPADDED genome length)"

echo ""
echo "============================================================"
echo "DONE. Files written to this folder:"
echo "   1. $REPORT   (diagnostic report)"
echo "   2. prelim_filtered_depth.txt   (unpadded control depth)"
echo ""
echo "To compare the control with the corrected profile, copy"
echo "prelim_filtered_depth.txt next to 02_fold_and_plot.py (with"
echo "padded_depth.txt and padding_metadata.txt) and run it."
echo "============================================================"
