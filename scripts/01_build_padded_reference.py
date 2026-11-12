#!/usr/bin/env python3
"""
01_build_padded_reference.py

Purpose
-------
Builds the "padded" linear reference used to correct for circularity.
Takes the circularized mitogenome FASTA and appends the first K bp of the
sequence onto the end, producing a reference of length N+K. This gives
read pairs that truly span the biological origin (the point where the
circle was cut to write a linear FASTA) somewhere contiguous to map to.

Also writes a small metadata file (padding_metadata.txt) recording the
original length (N) and the padding length (K), which the downstream
folding/plotting script needs to reconstruct true per-base depth. For
provenance, the metadata file also records the pipeline version and the MD5
checksum of the assembly FASTA, which 02_fold_and_plot.py copies into the
reported statistics.

Second mode: padding diagnostic report
--------------------------------------
After the padded mapping, run_pipeline.sh calls this script again with
--diagnose (Phase 6). Run inside a coverage_run/ folder, it reads
padding_metadata.txt, padded_depth.txt and filtered.bam and writes
diagnostic_report_<specimen>_<date>.txt: the padded geometry, alignment
record counts, a test for reads counted twice (expected 0), MAPQ
distributions in both copies of the duplicated window and the interior,
and soft-clipping near the junction. It needs samtools on the PATH. This
replaces the separate diagnose_padding_artifact.sh of earlier versions.

Usage
-----
    python3 01_build_padded_reference.py assembly.fasta K padded_assembly.fasta metadata_out.txt
    python3 01_build_padded_reference.py --diagnose [coverage_run_folder]
"""

import hashlib
import os
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime

__version__ = "4.5"   # coverage pipeline version

# Valid IUPAC nucleotide codes (upper-case; sequence is upper-cased on read).
VALID_BASES = set("ACGTUNRYSWKMBDHV")


def read_single_fasta(path):
    header = None
    seq_chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line[1:]
                else:
                    break  # only use first record
            else:
                seq_chunks.append(line.strip())
    if header is None:
        raise ValueError(f"No FASTA records found in {path}")
    return header, "".join(seq_chunks).upper()


def validate_sequence(seq, path):
    """Check every character is a valid nucleotide/IUPAC ambiguity code.
    Non-nucleotide characters get counted toward the sequence length here
    but are not accepted by bowtie2-build/bwa, silently producing a shorter
    indexed reference than this script assumes -- which shows up downstream
    as a confusing depth-file row-count mismatch. Catching it here,
    immediately, is much clearer.

    The most common case in a NOVOPlasty assembly is '*'. Per the NOVOPlasty
    documentation, '*' is a deliberate flag meaning the nucleotide BEFORE it
    is a possible insertion/deletion -- typically where the length of a
    single-nucleotide repeat could not be determined exactly because of
    systemic Illumina errors or a repetitive region. It is NOT a stray
    character and NOT equivalent to 'N'.
    (https://github.com/ndierckx/NOVOPlasty/wiki/Interpretation-&-post-processing)
    """
    bad_positions = [i + 1 for i, ch in enumerate(seq) if ch not in VALID_BASES]
    if bad_positions:
        shown = bad_positions[:20]
        more = f" (+{len(bad_positions) - 20} more)" if len(bad_positions) > 20 else ""
        bad_chars = sorted(set(seq[p - 1] for p in bad_positions))
        msg = (
            f"{path} contains {len(bad_positions)} character(s) that are not "
            f"nucleotides/IUPAC ambiguity codes: {bad_chars}. "
            f"Position(s) (1-based): {shown}{more}."
        )
        if "*" in bad_chars:
            msg += (
                "\n\n'*' is a NOVOPlasty uncertainty flag: the nucleotide "
                "immediately BEFORE each '*' is a possible insertion or "
                "deletion (single-nucleotide-repeat length could not be "
                "resolved). NOVOPlasty advises resolving these manually or "
                "deleting them before downstream analysis.\n"
                "  - Resolving manually (preferred): inspect the read pileup "
                "at these positions and set the repeat to its true length. "
                "Important if the site falls in a coding region, since an "
                "unresolved indel shifts the reading frame.\n"
                "  - Deleting: acceptable per the NOVOPlasty docs, but note "
                "each deletion shifts all downstream coordinates by one and "
                "commits to the assumption that the flagged base was spurious.\n"
                "  - Do NOT substitute 'N': '*' flags the preceding base, it "
                "is not itself a base placeholder."
            )
        else:
            msg += (
                "\n\nRemove or resolve these before re-running; they are not "
                "accepted by bowtie2-build and will silently shorten the "
                "indexed reference."
            )
        raise ValueError(msg)


def md5_of_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_fasta(path, header, seq, width=70):
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")


# ============================ DIAGNOSTIC REPORT ==============================

def _samtools(*args):
    return subprocess.run(["samtools", *args], check=True, capture_output=True,
                          text=True, errors="replace").stdout


def diagnose(folder="."):
    """Write the padding diagnostic report for a coverage_run/ folder.
    Returns 0 on success, 1 if inputs or samtools are missing."""
    folder = os.path.abspath(folder)
    missing = [f for f in ("padding_metadata.txt", "padded_depth.txt", "filtered.bam")
               if not os.path.isfile(os.path.join(folder, f))]
    if missing:
        print(f"MISSING in {folder}: {', '.join(missing)}")
        print("Run this inside a specimen's coverage_run/ folder. Nothing was written.")
        return 1
    if shutil.which("samtools") is None:
        print("ERROR: samtools not found on PATH. Run this inside WSL, not Windows.")
        return 1

    meta = {}
    with open(os.path.join(folder, "padding_metadata.txt")) as fh:
        meta_text = fh.read()
    for line in meta_text.splitlines():
        parts = line.split("\t")
        if len(parts) == 2:
            meta[parts[0]] = parts[1]
    k = int(meta["padding_length_K"])
    with open(os.path.join(folder, "padded_depth.txt")) as fh:
        total = sum(1 for _ in fh)
    n = total - k
    bam = os.path.join(folder, "filtered.bam")

    # One pass over the alignments collects everything below.
    window_copy = {}                 # read (QNAME + mate) -> set of copies
    mapq_a, mapq_b, mapq_int = Counter(), Counter(), Counter()
    n_int = 0
    cigars = Counter()
    proc = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE,
                            text=True, errors="replace")
    for line in proc.stdout:
        f = line.split("\t", 6)
        flag, pos, mapq, cigar = int(f[1]), int(f[3]), f[4], f[5]
        mate = "/1" if flag & 64 else ("/2" if flag & 128 else "/0")
        if pos <= k:
            window_copy.setdefault(f[0] + mate, set()).add("A")
            mapq_a[int(mapq)] += 1
        elif pos > n:
            window_copy.setdefault(f[0] + mate, set()).add("B")
            mapq_b[int(mapq)] += 1
        elif n_int < 200000:
            mapq_int[int(mapq)] += 1
            n_int += 1
        if pos > n - 200:
            cigars[cigar] += 1
    proc.wait()
    dup_reads = sum(1 for v in window_copy.values() if len(v) == 2)

    # Name the report after the specimen folder: the parent of coverage_run/
    # when run by the driver, otherwise the folder itself (by-hand runs).
    if os.path.basename(folder) == "coverage_run":
        specimen = os.path.basename(os.path.dirname(folder))
    else:
        specimen = os.path.basename(folder)
    report = os.path.join(folder, f"diagnostic_report_{specimen}_{datetime.now():%Y%m%d}.txt")
    sq = [l for l in _samtools("view", "-H", bam).splitlines() if l.startswith("@SQ")]

    def table(counter):
        return [f"    {c:>8}  {q}" for q, c in sorted(counter.items())] or ["    (none)"]

    out = [
        "=" * 60,
        f"PADDING DIAGNOSTIC (pipeline version {__version__})",
        f"Specimen folder : {specimen}",
        f"Path            : {folder}",
        f"Date            : {datetime.now():%Y-%m-%d %H:%M:%S}",
        f"samtools        : {_samtools('--version').splitlines()[0]}",
        "=" * 60, "",
        "--- geometry ---",
        f"padded_depth.txt rows (N+K) : {total}",
        f"K (padding length)          : {k}",
        f"N (true genome length)      : {n}",
        "", "padding_metadata.txt contents:",
        *["    " + l for l in meta_text.splitlines()],
        "", f"--- BAM header reference length (should equal N+K = {total}) ---",
        *["    " + l for l in sq],
        "", "--- alignment record counts in filtered.bam ---",
        f"total records          : {_samtools('view', '-c', bam).strip()}",
        f"secondary (flag 256)   : {_samtools('view', '-c', '-f', '256', bam).strip()}",
        f"supplementary (2048)   : {_samtools('view', '-c', '-f', '2048', bam).strip()}",
        f"duplicates (1024)      : {_samtools('view', '-c', '-f', '1024', bam).strip()}",
        "", "--- DOUBLE-COUNT TEST ---",
        "Individual READS (QNAME + mate number) with alignments in BOTH copies",
        f"of the duplicated window (copy 1 = pos 1..{k}, copy 2 = pos > {n}).",
        "Mate number is included because both mates of a pair share a QNAME.",
        "If secondary/supplementary counts above are 0, this MUST be 0: with one",
        "alignment record per read, a read cannot contribute depth twice.",
        f"reads present in BOTH copies : {dup_reads}",
        "", f"--- MAPQ distribution: window copy 1 (pos 1..{k}) ---", "       count  MAPQ",
        *table(mapq_a),
        "", f"--- MAPQ distribution: window copy 2 (pos > {n}) ---", "       count  MAPQ",
        *table(mapq_b),
        "", f"--- MAPQ distribution: interior (pos {k + 1}..{n}, first 200k records) ---",
        "       count  MAPQ",
        *table(mapq_int),
        "", f"--- CIGAR soft-clipping near the junction (pos > {n - 200}) ---",
        "       count  CIGAR (top 15)",
        *[f"    {c:>8}  {cg}" for cg, c in cigars.most_common(15)],
        "", "=" * 60, "END OF DIAGNOSTIC", "=" * 60,
    ]
    with open(report, "w") as fh:
        fh.write("\n".join(out) + "\n")
    print("\n".join(out))
    print(f"\nDiagnostic report written to {report}")
    return 0


def main():
    if len(sys.argv) >= 2 and sys.argv[1] == "--diagnose":
        sys.exit(diagnose(sys.argv[2] if len(sys.argv) > 2 else "."))
    if len(sys.argv) != 5:
        print("Usage: python3 01_build_padded_reference.py "
              "assembly.fasta K padded_assembly.fasta metadata_out.txt\n"
              "       python3 01_build_padded_reference.py --diagnose [coverage_run_folder]",
              file=sys.stderr)
        sys.exit(1)

    assembly_path, k_str, padded_out_path, metadata_out_path = sys.argv[1:5]
    k = int(k_str)

    header, seq = read_single_fasta(assembly_path)
    validate_sequence(seq, assembly_path)
    n = len(seq)

    if k <= 0:
        raise ValueError("K (padding length) must be a positive integer")
    if k >= n:
        raise ValueError(
            f"K ({k}) must be smaller than the genome length ({n}). "
            f"Something is off with the insert-size estimate -- check "
            f"prelim.stats.txt before continuing."
        )

    padded_seq = seq + seq[:k]
    padded_header = f"{header}|padded_K{k}"

    write_fasta(padded_out_path, padded_header, padded_seq)

    with open(metadata_out_path, "w") as fh:
        fh.write(f"original_header\t{header}\n")
        fh.write(f"original_length_N\t{n}\n")
        fh.write(f"padding_length_K\t{k}\n")
        fh.write(f"padded_length\t{n + k}\n")
        fh.write(f"assembly_file\t{assembly_path}\n")
        fh.write(f"assembly_md5\t{md5_of_file(assembly_path)}\n")
        fh.write(f"pipeline_version\t{__version__}\n")

    print(f"Original length (N): {n} bp")
    print(f"Padding length (K):  {k} bp")
    print(f"Padded reference written to: {padded_out_path} ({n + k} bp)")
    print(f"Metadata written to: {metadata_out_path}")


if __name__ == "__main__":
    main()
