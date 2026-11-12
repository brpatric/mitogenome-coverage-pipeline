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

Usage
-----
    python3 01_build_padded_reference.py assembly.fasta K padded_assembly.fasta metadata_out.txt
"""

import hashlib
import sys

__version__ = "4.2"   # coverage pipeline version

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


def main():
    if len(sys.argv) != 5:
        print("Usage: python3 01_build_padded_reference.py "
              "assembly.fasta K padded_assembly.fasta metadata_out.txt",
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
