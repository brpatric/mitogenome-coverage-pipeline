#!/usr/bin/env python3
"""
00_check_terminal_overlap.py

Purpose
-------
Before padding a circularized mitogenome assembly for coverage mapping, we
need to confirm the FASTA does NOT already contain a duplicated stretch at
the start/end (some circularizers leave the overlap that was used to close
the loop instead of trimming it). If that redundancy is present and we pad
on top of it, coverage in that region will be double-counted.

This script checks for an exact suffix == prefix match and reports the
longest one found at or above a minimum length, if any. It does not modify
your file -- it only reports, so a human can decide how to handle a hit.

Why there is a minimum length
-----------------------------
A short suffix/prefix match means nothing. The last base of a sequence
matches the first base by chance about one time in four; a 2 bp match about
one time in sixteen, and so on. Only a match long enough to be implausible
by chance indicates genuine leftover circularization overlap, which in
practice is tens to hundreds of bases. The default minimum is 20 bp, at
which a chance match is about one in a trillion for random sequence.
Matches shorter than the minimum are reported as an informational note and
do NOT count as redundancy.

Usage
-----
    python3 00_check_terminal_overlap.py assembly.fasta [--window 2000]
                                                        [--min-overlap 20]

Exit codes
----------
    0 = no redundancy at or above the minimum length (safe to proceed)
    1 = redundancy found at or above the minimum length (investigate)
"""

import argparse
import sys

__version__ = "4.1.1"   # coverage pipeline version

# Valid IUPAC nucleotide codes (upper-case; sequence is upper-cased on read).
VALID_BASES = set("ACGTUNRYSWKMBDHV")


def read_single_fasta(path):
    """Minimal FASTA reader. Returns (header, sequence) for the first record.
    Warns if more than one record is present (only the first is used)."""
    header = None
    seq_chunks = []
    n_records = 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                n_records += 1
                if n_records == 1:
                    header = line[1:]
                elif n_records == 2:
                    # stop reading once we hit a second record
                    break
            else:
                if n_records == 1:
                    seq_chunks.append(line.strip())
    if header is None:
        raise ValueError(f"No FASTA records found in {path}")
    if n_records > 1:
        print(f"WARNING: {path} contains {n_records}+ records. "
              f"Only the first record ('{header}') was checked. "
              f"Confirm this is your circularized mitogenome contig.",
              file=sys.stderr)
    return header, "".join(seq_chunks).upper()


def longest_suffix_prefix_overlap(seq, window, min_overlap=1):
    """Return the length of the longest exact match between seq[-k:] and
    seq[:k], searching k from `window` down to `min_overlap`. Returns 0 if
    no match of at least `min_overlap` bases is found."""
    window = min(window, len(seq) // 2)
    if min_overlap < 1:
        min_overlap = 1
    for k in range(window, min_overlap - 1, -1):
        if seq[-k:] == seq[:k]:
            return k
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fasta", help="Path to the circularized assembly FASTA")
    ap.add_argument("--window", type=int, default=2000,
                     help="Max bp to check at each end (default: 2000)")
    ap.add_argument("--min-overlap", type=int, default=20, dest="min_overlap",
                     help="Shortest match treated as real redundancy, in bp "
                          "(default: 20). Shorter matches occur by chance and "
                          "are reported as a note only.")
    args = ap.parse_args()

    header, seq = read_single_fasta(args.fasta)
    n = len(seq)
    print(f"00_check_terminal_overlap.py (pipeline version {__version__})")
    print(f"Sequence: {header}")
    print(f"Length: {n} bp")

    bad_positions = [i + 1 for i, ch in enumerate(seq) if ch not in VALID_BASES]
    if bad_positions:
        shown = bad_positions[:20]
        more = f" (+{len(bad_positions) - 20} more)" if len(bad_positions) > 20 else ""
        bad_chars = sorted(set(seq[p - 1] for p in bad_positions))
        print("")
        print(f"WARNING: {len(bad_positions)} non-nucleotide character(s) found: {bad_chars}")
        print(f"Position(s) (1-based): {shown}{more}")
        if "*" in bad_chars:
            print("'*' is a NOVOPlasty uncertainty flag: the nucleotide immediately")
            print("BEFORE each '*' is a possible insertion/deletion (single-nucleotide-")
            print("repeat length could not be resolved). It is NOT a stray character and")
            print("NOT equivalent to 'N'. NOVOPlasty advises resolving these manually")
            print("(inspect the read pileup and set the true repeat length) or deleting")
            print("them. Resolve before finalizing -- an unresolved indel inside a coding")
            print("region shifts the reading frame.")
        print("Either way they must not remain in the FASTA used for mapping: they are")
        print("rejected by bowtie2-build, silently shrinking the indexed reference and")
        print("causing a depth-file length mismatch later.")
        print("")

    print(f"Checking for exact suffix/prefix redundancy between "
          f"{args.min_overlap} and {min(args.window, n // 2)} bp...")

    overlap = longest_suffix_prefix_overlap(seq, args.window, args.min_overlap)

    if overlap == 0:
        # Nothing at or above the threshold. Report any shorter match, but
        # describe it honestly: a 1-5 bp match really is chance, while a
        # 15 bp match is not, even though both are below the default cutoff.
        trivial = longest_suffix_prefix_overlap(seq, args.min_overlap - 1, 1)
        print(f"RESULT: No terminal redundancy of {args.min_overlap} bp or more "
              f"detected.")
        if 0 < trivial <= 5:
            print(f"NOTE: the last {trivial} bp match the first {trivial} bp. A "
                  f"match this short is expected by chance (about 1 in "
                  f"{4 ** trivial} for random sequence) and is NOT evidence of "
                  f"leftover circularization overlap. No action needed.")
        elif trivial > 5:
            print(f"NOTE: the last {trivial} bp match the first {trivial} bp. "
                  f"This is below the {args.min_overlap} bp reporting threshold, "
                  f"but it is longer than chance alone would readily explain "
                  f"(about 1 in {4 ** trivial:,} for random sequence). It may be "
                  f"a short real overlap or a genuine repeat at the origin. "
                  f"Worth a look before you proceed; re-run with "
                  f"--min-overlap {trivial} to treat it as redundancy.")
        print("You can proceed directly to padding this assembly.")
    else:
        print(f"RESULT: Exact terminal redundancy of {overlap} bp detected.")
        print(f"The last {overlap} bp of the sequence exactly match the first {overlap} bp.")
        print("A match this long is very unlikely by chance, so this is probably")
        print("leftover overlap from circularization that was not trimmed.")
        print("Recommended fix: remove the last "
              f"{overlap} bp from the sequence before padding, so the two ends")
        print("represent a single, non-redundant copy of the junction.")
        print("")
        print("If you confirm this is real redundancy (not a biological repeat,")
        print("which can also produce a long terminal match), trim it and re-run")
        print("this check on the trimmed FASTA before continuing to the padding step.")

    return 0 if overlap == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
