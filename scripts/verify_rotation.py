#!/usr/bin/env python3
"""
verify_rotation.py

Purpose
-------
Confirms that a manually rotated circular assembly (e.g. one you rotated so
that COI comes first) is an EXACT rotation of the original assembler output
-- no bases duplicated, lost, or altered at the seam where you cut and
rejoined the circle.

This is a stronger, more direct test than the terminal-overlap heuristic.
A genuine rotation of a circular sequence must satisfy all three:

    1. identical length
    2. identical base composition
    3. the rotated sequence appears within (original + original)

Condition 3 is the decisive one: it holds if and only if the rotated
sequence is some cyclic shift of the original. Duplicating or dropping even
a single base at the seam breaks it.

Usage
-----
    python3 verify_rotation.py original_from_assembler.fasta rotated_COIfirst.fasta

Exit code 0 = verified clean rotation, 1 = problem found.
"""

import sys
from collections import Counter

__version__ = "4.2.2"   # coverage pipeline version


def read_single_fasta(path):
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line[1:]
                else:
                    break
            else:
                chunks.append(line.strip())
    if header is None:
        raise ValueError(f"No FASTA records found in {path}")
    return header, "".join(chunks).upper()


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 verify_rotation.py original.fasta rotated.fasta",
              file=sys.stderr)
        sys.exit(2)

    orig_path, rot_path = sys.argv[1], sys.argv[2]
    orig_hdr, orig = read_single_fasta(orig_path)
    rot_hdr, rot = read_single_fasta(rot_path)

    print(f"verify_rotation.py (pipeline version {__version__})")
    print(f"original : {orig_hdr}")
    print(f"           {orig_path}  ({len(orig)} bp)")
    print(f"rotated  : {rot_hdr}")
    print(f"           {rot_path}  ({len(rot)} bp)")
    print("")

    ok = True

    # --- 1. length ---
    if len(orig) == len(rot):
        print(f"[PASS] length identical ({len(orig)} bp)")
    else:
        diff = len(rot) - len(orig)
        print(f"[FAIL] length differs by {diff:+d} bp "
              f"(original {len(orig)}, rotated {len(rot)})")
        print(f"       {'Extra bases' if diff > 0 else 'Missing bases'} "
              f"at the seam is the usual cause of a {abs(diff)} bp difference.")
        ok = False

    # --- 2. composition ---
    co, cr = Counter(orig), Counter(rot)
    if co == cr:
        print("[PASS] base composition identical")
    else:
        print("[FAIL] base composition differs:")
        for base in sorted(set(co) | set(cr)):
            if co.get(base, 0) != cr.get(base, 0):
                print(f"       {base}: original={co.get(base,0)}  rotated={cr.get(base,0)}")
        ok = False

    # --- 3. cyclic-shift test (decisive) ---
    if len(orig) == len(rot) and rot in (orig + orig):
        offset = (orig + orig).index(rot)
        print(f"[PASS] rotated sequence is an exact cyclic shift of the original")
        print(f"       rotation offset: {offset} bp "
              f"(rotated position 1 == original position {offset + 1})")
    else:
        print("[FAIL] rotated sequence is NOT an exact cyclic shift of the original")
        print("       The sequence content changed at the seam -- this is not a")
        print("       pure rotation. Re-do the rotation from the assembler output.")
        ok = False

    print("")
    if ok:
        print("VERDICT: clean rotation. Safe to use for coverage mapping.")
    else:
        print("VERDICT: PROBLEM FOUND. Do not use this file for coverage mapping")
        print("         until the rotation is redone.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
