#!/usr/bin/env python3
"""
02_fold_and_plot.py

Author
------
L. Brian Patrick, Department of Biological Sciences, Dakota Wesleyan
University, Mitchell, SD, USA. ORCID 0000-0003-1969-9600 (https://orcid.org/0000-0003-1969-9600)
Repository: https://github.com/brpatric/mitogenome-coverage-pipeline

Purpose
-------
Takes the raw per-base depth from the PADDED reference (length N+K) and
"folds" it back onto the true circular genome (length N):

    folded[pos] = padded_depth[pos]                      for pos in (K, N]
    folded[pos] = padded_depth[pos] + padded_depth[N+pos] for pos in [1, K]

This adds back, onto positions 1..K, the depth from reads/pairs that were
placed in the duplicated tail because they truly span the biological
origin. The result is a per-base depth profile of the real genome, which is
what gets reported and plotted.

Two checks follow, each run when its input file is present:

  1. Interior check (prelim_filtered_depth.txt). Compares the folded
     profile with the unpadded control across the genome interior --
     positions K+1 .. N-K, which padding does not touch and where the
     control is free of the origin artifact.

  2. End-window check (origin_check_depth.txt, pipeline version 4.2+).
     The interior check cannot see the ends, which is exactly where the
     padding acts. Phase 8 of run_pipeline.sh therefore maps the reads to a
     copy of the assembly cut at the opposite side of the circle, where the
     original start/end region lies deep in the interior. This script
     compares that profile with the folded profile across the end windows
     (the last W and first W positions, W = K or less) in 200-bp bins.
     If no bin differs by 10% or more AND 3x or more, the check passes and
     the folded profile is reported. If any bin does, the check fails: the
     end windows of the reported profile are taken from the origin-check
     mapping instead, the reason is written to the outputs, and the folded
     profile is kept separately for the record.

If either input file is absent, that check is skipped with a note and
everything else runs as normal.

How to run
----------
Two ways:

  1. Automatically, by run_pipeline.sh (Phase 9, version 4.4+). If this
     script sits in the same folder as run_pipeline.sh, the driver runs it
     at the end of every run, reading the files in coverage_run/ and
     writing the outputs to coverage_results/ in the folder the driver was
     run from. The prefix and plot title come from the driver's settings.

  2. By hand, from the Windows Python IDLE shell (Run > Run Module, or F5)
     or from WSL, with no command-line arguments. Put these files in the
     SAME FOLDER as this script (all from coverage_run/):
    padded_depth.txt            required
    padding_metadata.txt        required
    prelim_filtered_depth.txt   optional, needed for the interior check
    origin_check_depth.txt      optional, needed for the end-window check
    run_provenance.txt          optional, adds the mapping details to the report
When you run it, the program will ask you to type an output prefix (e.g.
MySpecies_sample01) -- this is used to name all the output files.

It only needs pandas, matplotlib, and numpy -- no bioinformatics tools --
so it also runs fine inside WSL if you'd rather use that instead.

Outputs (by hand: next to this script; from the driver: coverage_results/)
-----------------------------------------------------------------------
    <prefix>_folded_depth.txt        per-base depth table reported (Pos, Depth)
    <prefix>_coverage_stats.txt      mean/sd/min/max/breadth of that profile
    <prefix>_linear_coverage.png     linear coverage plot
    <prefix>_circular_coverage.png   circular/polar coverage plot
    <prefix>_control_comparison.txt  interior check
                                     (only if prelim_filtered_depth.txt is present)
    <prefix>_end_window_check.txt    end-window check
                                     (only if origin_check_depth.txt is present)
    <prefix>_folded_depth_before_end_window_correction.txt
                                     the folded profile, kept for the record
                                     (only if the end-window check failed)
    <prefix>_report.txt              everything to report in the Methods and
                                     Results, with draft sentences. Mapping
                                     details (software versions, insert size,
                                     rotation and overlap checks) come from
                                     run_provenance.txt if it is with the inputs

03_compare_to_control.py produces the same interior comparison as a
standalone step, for depth profiles folded with an earlier version of this
script.

Provenance: the statistics and check files end with a provenance block
recording the pipeline version, date, Python and library versions, and MD5
checksums of the input files (and, via padding_metadata.txt, of the
assembly). The values above the block are unchanged by it.
"""

import argparse
import hashlib
import os
import sys
import math
import platform
from datetime import datetime

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

__version__ = "4.5"   # coverage pipeline version

# Archived releases (Zenodo), for the citation in <prefix>_report.txt:
# version -> (year, DOI). Add each new release here once it has a DOI.
PIPELINE_RELEASES = {"4.5": ("2026", "10.5281/zenodo.24181781")}
REPO_URL = "https://github.com/brpatric/mitogenome-coverage-pipeline"

# References cited in the draft text of <prefix>_report.txt.
REFERENCES = [
    "Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, "
    "Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. 2021. Twelve years of "
    "SAMtools and BCFtools. GigaScience 10(2): giab008. "
    "doi:10.1093/gigascience/giab008",
    "Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, "
    "Cournapeau D, Wieser E, Taylor J, Berg S, Smith NJ, Kern R, Picus M, "
    "Hoyer S, van Kerkwijk MH, Brett M, Haldane A, Fernandez del Rio J, Wiebe M, "
    "Peterson P, Gerard-Marchant P, Sheppard K, Reddy T, Weckesser W, Abbasi H, "
    "Gohlke C, Oliphant TE. 2020. Array programming with NumPy. Nature "
    "585(7825): 357-362. doi:10.1038/s41586-020-2649-2",
    "Hunter JD. 2007. Matplotlib: a 2D graphics environment. Computing in "
    "Science & Engineering 9(3): 90-95. doi:10.1109/MCSE.2007.55",
    "Langmead B, Salzberg SL. 2012. Fast gapped-read alignment with Bowtie 2. "
    "Nature Methods 9(4): 357-359. doi:10.1038/nmeth.1923",
    "McKinney W. 2010. Data structures for statistical computing in Python. "
    "In: van der Walt S, Millman J, editors. Proceedings of the 9th Python in "
    "Science Conference. p. 56-61. doi:10.25080/Majora-92bf1922-00a",
]
ASSEMBLER_REFERENCES = [
    "Dierckxsens N, Mardulyn P, Smits G. 2017. NOVOPlasty: de novo assembly of "
    "organelle genomes from whole genome data. Nucleic Acids Research 45(4): "
    "e18. doi:10.1093/nar/gkw955",
    "Jin JJ, Yu WB, Yang JB, Song Y, dePamphilis CW, Yi TS, Li DZ. 2020. "
    "GetOrganelle: a fast and versatile toolkit for accurate de novo assembly "
    "of organelle genomes. Genome Biology 21(1): 241. "
    "doi:10.1186/s13059-020-02154-5",
]


def pipeline_citation(version):
    """In-text citation and reference-list entry for this pipeline version."""
    year, doi = PIPELINE_RELEASES.get(version, (None, None))
    if doi:
        return (f"Patrick {year}",
                f"Patrick LB. {year}. Mitogenome Coverage Pipeline, version "
                f"{version}. Zenodo. doi:{doi}. {REPO_URL}")
    return ("Patrick [year]",
            f"Patrick LB. [year]. Mitogenome Coverage Pipeline, version {version}. "
            f"Zenodo. doi:[DOI of the version {version} release; see {REPO_URL}]. "
            f"{REPO_URL}")

# ============================= EDIT THESE ===================================
PLOT_TITLE = "Per-base Read Coverage"     # set to the full title you want shown, e.g. "Per-base Read Coverage: Genus species"
LINEAR_COLOR = "#648FFF"
MEAN_LINE_COLOR = "#FE6100"
# =============================================================================

# End-window check settings (rarely need changing). A bin fails only if it
# differs by at least END_MAX_PCT percent AND at least END_MIN_ABS reads.
END_BIN_BP = 200          # bin width for the comparison, bp
END_MAX_PCT = 10.0        # percent difference that fails a bin
END_MIN_ABS = 3.0         # absolute difference (x) that fails a bin
END_MARGIN_BP = 500       # keep the window this far from the origin-check
                          # mapping's own sequence ends

# By hand, inputs and outputs are in the same folder as this script. The
# driver passes --input-dir and --output-dir instead.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
IN_DIR = SCRIPT_DIR
OUT_DIR = SCRIPT_DIR
DEPTH_FILENAME = "padded_depth.txt"
METADATA_FILENAME = "padding_metadata.txt"
CONTROL_FILENAME = "prelim_filtered_depth.txt"
ORIGIN_FILENAME = "origin_check_depth.txt"
PROVENANCE_FILENAME = "run_provenance.txt"   # optional: mapping details for the report


def read_metadata(path):
    meta = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                meta[parts[0]] = parts[1]
    return meta


def md5_of_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance_lines(script_name, input_paths, meta):
    """Provenance block appended to each reported output file: what
    produced it, when, with which library versions, and from which inputs
    (MD5 checksums). The values above the block are unaffected."""
    lines = [
        "",
        "--- Provenance ---",
        f"Produced by       : {script_name} (pipeline version {__version__})",
        f"Date              : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Python            : {platform.python_version()}",
        f"NumPy / pandas / Matplotlib : {np.__version__} / {pd.__version__} / {matplotlib.__version__}",
    ]
    if "assembly_md5" in meta:
        lines.append(f"Assembly (MD5)    : {meta.get('assembly_file', 'assembly')} {meta['assembly_md5']}")
    lines.append(f"Padding built by  : pipeline version {meta.get('pipeline_version', 'earlier than 4.1 (not recorded)')}")
    for p in input_paths:
        lines.append(f"Input (MD5)       : {os.path.basename(p)} {md5_of_file(p)}")
    return lines


def read_depth_column(path):
    return pd.read_csv(path, sep="\t", header=None,
                       names=["Ref", "Pos", "Depth"])["Depth"].to_numpy()


def fold_depth(padded_depth, n, k):
    """padded_depth: 1D numpy array, 0-indexed array where index i
    corresponds to 1-based genomic position i+1 on the padded reference.
    Returns folded depth array of length n (index i -> position i+1)."""
    if len(padded_depth) != n + k:
        raise ValueError(
            f"padded_depth has {len(padded_depth)} rows but metadata says "
            f"N+K = {n + k}. Did you point this at the right depth file?"
        )
    folded = padded_depth[:n].copy()
    # positions 1..K (0-indexed 0..K-1) get the padded tail added back
    folded[:k] = folded[:k] + padded_depth[n:n + k]
    return folded


def compare_to_control(folded, control, n, k, folded_label):
    """Compare corrected and control profiles across the interior
    (1-based positions K+1 .. N-K). Returns the report as a list of lines,
    or None (with a printed note) if the comparison cannot be made.
    The comparison lines are identical to those of 03_compare_to_control.py."""
    if len(control) != n:
        print(f"NOTE: {CONTROL_FILENAME} has {len(control)} positions, but "
              f"metadata says N = {n}. Files are from different runs? "
              f"Control comparison skipped.")
        return None
    if n - 2 * k < 1:
        print(f"NOTE: genome (N = {n}) too short for an interior with "
              f"K = {k}. Control comparison skipped.")
        return None

    # interior = 1-based positions K+1 .. N-K  ->  0-based slice [k : n-k]
    lo, hi = k + 1, n - k
    f_int = folded[k:n - k]
    c_int = control[k:n - k]

    f_mean = float(np.mean(f_int))
    c_mean = float(np.mean(c_int))
    pct_diff = 100 * (f_mean - c_mean) / c_mean

    f_min = int(np.min(f_int))
    c_min = int(np.min(c_int))
    f_min_pos = np.flatnonzero(f_int == f_min) + lo
    c_min_pos = np.flatnonzero(c_int == c_min) + lo
    same_min = (f_min == c_min) and np.array_equal(f_min_pos, c_min_pos)

    def fmt_pos(arr):
        shown = ", ".join(str(p) for p in arr[:10])
        return shown + (f" (+{len(arr) - 10} more)" if len(arr) > 10 else "")

    return [
        f"Folded profile : {folded_label}",
        f"Control profile: {CONTROL_FILENAME}",
        f"N = {n} bp, K = {k} bp",
        f"Interior compared: positions {lo}-{hi} ({hi - lo + 1} bp)",
        "",
        f"Interior mean depth, corrected : {f_mean:.4f}",
        f"Interior mean depth, control   : {c_mean:.4f}",
        f"Difference (corrected - control): {pct_diff:+.3f}%",
        "",
        f"Interior minimum, corrected : {f_min} at position(s) {fmt_pos(f_min_pos)}",
        f"Interior minimum, control   : {c_min} at position(s) {fmt_pos(c_min_pos)}",
        f"Same minimum at same position(s): {'YES' if same_min else 'NO'}",
    ]


def end_window_check(folded, origin_shifted, n, k, s):
    """Compare the folded profile with the origin-check profile across the
    end windows (circularly contiguous region centred on the sequence
    start/end). origin_shifted is the depth from the mapping to the
    assembly rotated so that it starts at original position s+1.

    Returns (lines, passed, corrected) or None if the check cannot run.
    corrected is the folded profile with the end windows replaced by the
    origin-check profile (used only if the check fails)."""
    if len(origin_shifted) != n:
        print(f"NOTE: {ORIGIN_FILENAME} has {len(origin_shifted)} positions, but "
              f"metadata says N = {n}. Files are from different runs? "
              f"End-window check skipped.")
        return None
    # Back to the assembly's own coordinates:
    # rotated index j holds original index (j + s) mod n.
    origin = np.roll(origin_shifted, s)

    # The origin-check mapping has its own sequence ends at original
    # position s / s+1; keep the window END_MARGIN_BP away from them.
    w = min(k, s - END_MARGIN_BP, (n - s) - END_MARGIN_BP)
    if w < END_BIN_BP:
        print(f"NOTE: genome (N = {n}) too short for an end-window check "
              f"with K = {k}. End-window check skipped.")
        return None

    idx = np.r_[np.arange(n - w, n), np.arange(0, w)]   # last w, then first w
    f_w = folded[idx].astype(float)
    o_w = origin[idx].astype(float)

    def pos_label(i):
        return idx[i] + 1

    failed_bins = []
    worst = None
    starts = list(range(0, len(idx), END_BIN_BP))
    if len(starts) > 1 and len(idx) - starts[-1] < END_BIN_BP // 2:
        starts.pop()          # fold a short last bin into the one before it
    ends = starts[1:] + [len(idx)]
    for b0, b1 in zip(starts, ends):
        fm, om = f_w[b0:b1].mean(), o_w[b0:b1].mean()
        diff = fm - om
        pct = 100 * diff / max(om, 1.0)
        entry = (pos_label(b0), pos_label(b1 - 1), fm, om, diff, pct)
        if worst is None or abs(pct) > abs(worst[5]):
            worst = entry
        if abs(diff) >= END_MIN_ABS and abs(pct) >= END_MAX_PCT:
            failed_bins.append(entry)

    passed = not failed_bins
    f_mean, o_mean = f_w.mean(), o_w.mean()
    per_base = np.abs(f_w - o_w)
    imax = int(np.argmax(per_base))
    step_f = abs(int(folded[0]) - int(folded[n - 1]))
    step_o = abs(int(origin[0]) - int(origin[n - 1]))

    corrected = folded.copy()
    corrected[idx] = origin[idx]

    lines = [
        f"Folded profile       : padding and folding (this run)",
        f"Origin-check profile : {ORIGIN_FILENAME} (assembly copy starting at "
        f"original position {s + 1})",
        f"N = {n} bp, K = {k} bp",
        f"End windows compared : positions {n - w + 1}-{n} and 1-{w} "
        f"({2 * w} bp, in {END_BIN_BP}-bp bins)",
        f"A bin fails if it differs by >= {END_MAX_PCT:g}% AND >= {END_MIN_ABS:g}x",
        "",
        f"End-window mean depth, folded       : {f_mean:.4f}",
        f"End-window mean depth, origin check : {o_mean:.4f}",
        f"Difference (folded - origin check)  : {100 * (f_mean - o_mean) / max(o_mean, 1):+.3f}%",
        f"Largest one-base difference         : {per_base[imax]:.0f}x at position {pos_label(imax)}",
        f"Step between position {n} and position 1: folded {step_f}x, origin check {step_o}x",
        f"Largest bin difference              : positions {worst[0]}-{worst[1]}: "
        f"folded {worst[2]:.1f}x, origin check {worst[3]:.1f}x ({worst[5]:+.1f}%)",
        f"Bins failing                        : {len(failed_bins)}",
    ]
    for e in failed_bins[:20]:
        lines.append(f"    positions {e[0]}-{e[1]}: folded {e[2]:.1f}x, "
                     f"origin check {e[3]:.1f}x ({e[4]:+.1f}x, {e[5]:+.1f}%)")
    if len(failed_bins) > 20:
        lines.append(f"    (+{len(failed_bins) - 20} more)")
    lines.append("")
    if passed:
        lines.append("RESULT: PASSED. The folded profile agrees with the origin-check")
        lines.append("profile across the end windows; the folded profile is reported.")
    else:
        lines.append("RESULT: FAILED. The folded profile departs from the origin-check")
        lines.append("profile near the sequence start/end. In the reported profile,")
        lines.append(f"positions {n - w + 1}-{n} and 1-{w} are taken from the origin-check")
        lines.append("mapping, where this region lies far from any sequence end. The")
        lines.append("folded profile is kept as")
        lines.append("<prefix>_folded_depth_before_end_window_correction.txt.")
    info = {"w": w, "passed": passed, "n_failed": len(failed_bins),
            "worst": worst, "f_mean": f_mean, "o_mean": o_mean,
            "step_f": step_f, "step_o": step_o, "s": s}
    return lines, passed, corrected, info


def read_provenance(path):
    """Key/value lines of run_provenance.txt ("Key : value")."""
    prov = {}
    if not os.path.isfile(path):
        return prov
    with open(path) as fh:
        for line in fh:
            if " : " in line:
                key, val = line.split(" : ", 1)
                prov.setdefault(key.strip(), val.strip())
    return prov


def build_report(prefix, profile, folded, control, n, k, meta, prov, end_info,
                 profile_note):
    """Everything to report, in one place: the Methods values, the Results
    values, and draft sentences built from them."""
    P = lambda key, default="(not recorded)": prov.get(key, default)
    have_prov = bool(prov)
    mean, sd = float(np.mean(profile)), float(np.std(profile))
    median = float(np.median(profile))
    dmin, dmax = int(np.min(profile)), int(np.max(profile))
    breadth = int(np.sum(profile > 0))
    pos_min = np.flatnonzero(profile == dmin) + 1
    pos_max = np.flatnonzero(profile == dmax) + 1

    def ranges(pos):
        """1-based positions as compact ranges, e.g. 4495-4505, 7001."""
        if len(pos) == 0:
            return "-"
        out, start, prev = [], pos[0], pos[0]
        for p in list(pos[1:]) + [None]:
            if p is not None and p == prev + 1:
                prev = p
                continue
            out.append(f"{start}" if start == prev else f"{start}-{prev}")
            if p is not None:
                start = prev = p
        shown = ", ".join(out[:5])
        return shown + (f" (+{len(out) - 5} more)" if len(out) > 5 else "")

    rotated = P("ROTATED", "")
    L = ["=" * 72,
         f"COVERAGE REPORT: {prefix}",
         f"Everything to report in the Methods and Results, from this run.",
         f"Produced by 02_fold_and_plot.py (pipeline version {__version__}), "
         f"{datetime.now():%Y-%m-%d %H:%M}",
         "=" * 72, ""]
    if not have_prov:
        L += ["NOTE: run_provenance.txt was not found with the input files, so the",
              "mapping details (software versions, insert size, rotation and overlap",
              "checks, alignment counts) are missing below. To include them, copy",
              "coverage_run/run_provenance.txt next to the other input files and run",
              "this again. (A one-command run with run_pipeline.sh includes them.)", ""]

    version = P('Pipeline version', __version__)
    cite_in, cite_ref = pipeline_citation(version)
    L += ["METHODS (procedure and software)",
          "-" * 72,
          f"Pipeline                    : Mitogenome Coverage Pipeline, version "
          f"{version}",
          f"Pipeline citation           : {cite_ref}",
          f"Assembly                    : {meta.get('original_header', '?')} ({n} bp)",
          f"Rotated after assembly      : {rotated or '(not recorded)'}",
          f"Overlap check, unrotated    : {P('Overlap (unrotated)', 'not run (assembly not rotated)' if rotated == 'no' else '(not recorded)')}",
          f"Rotation check              : {P('Rotation check')}",
          f"Overlap check, as mapped    : {P('Overlap check')}",
          f"Bowtie2                     : {P('bowtie2')} (--local mode)",
          f"SAMtools                    : {P('samtools')}",
          f"Alignment filter            : properly paired only (samtools view -f 2), no MAPQ threshold",
          f"Insert size (mean, SD)      : {P('Insert size mean')} bp, {P('Insert size SD')} bp",
          f"Padding length K            : {k} bp (mean insert + 4 SD, minimum 300 bp)",
          f"Python                      : {platform.python_version()}",
          f"NumPy / pandas / Matplotlib : {np.__version__} / {pd.__version__} / {matplotlib.__version__}",
          f"End-window check mapping    : assembly rotated by S = {meta.get('origin_check_shift_S', '?')} bp (half its length)",
          ""]

    L += ["RESULTS (values)",
          "-" * 72,
          f"Genome length N             : {n} bp",
          f"Padding length K            : {k} bp",
          f"Breadth of coverage         : {breadth} bp ({100 * breadth / n:.2f}%)",
          f"Mean depth +/- SD           : {mean:.1f} +/- {sd:.1f}x   (unrounded {mean:.4f} +/- {sd:.4f})",
          f"Median depth                : {median:.1f}x",
          f"Minimum depth               : {dmin}x at {ranges(pos_min)}",
          f"Maximum depth               : {dmax}x at {ranges(pos_max)}",
          f"Properly paired alignments  : {P('Properly paired alignments (padded)')} (padded mapping)",
          f"Reads counted twice         : {P('Diagnostic report').split('window: ')[-1].rstrip(')') if 'window: ' in P('Diagnostic report') else '(not recorded)'} (padding diagnostic report; expected 0)",
          ""]

    interior = None
    if control is not None and len(control) == n and n - 2 * k >= 1:
        f_int, c_int = folded[k:n - k], control[k:n - k]
        fm, cm = float(np.mean(f_int)), float(np.mean(c_int))
        pct = 100 * (fm - cm) / cm
        fmin, cmin = int(np.min(f_int)), int(np.min(c_int))
        fpos = np.flatnonzero(f_int == fmin) + k + 1
        cpos = np.flatnonzero(c_int == cmin) + k + 1
        same = fmin == cmin and np.array_equal(fpos, cpos)
        interior = (pct, same)
        L += ["Interior check (validation)",
              f"  Positions compared        : {k + 1}-{n - k} ({n - 2 * k} bp)",
              f"  Mean depth, corrected     : {fm:.4f}x",
              f"  Mean depth, unpadded      : {cm:.4f}x",
              f"  Difference                : {pct:+.3f}%",
              f"  Minimum, corrected        : {fmin}x at {ranges(fpos)}",
              f"  Minimum, unpadded         : {cmin}x at {ranges(cpos)}",
              f"  Same minimum, same place  : {'YES' if same else 'NO'}",
              "  (Interior means are NOT the specimen's mean depth; report the",
              "   whole-genome values above.)", ""]
    else:
        L += ["Interior check (validation)  : not run", ""]

    if end_info is not None:
        w, worst = end_info["w"], end_info["worst"]
        L += ["End-window check (validation)",
              f"  Windows compared          : positions {n - w + 1}-{n} and 1-{w} (200-bp bins)",
              f"  Result                    : {'PASSED' if end_info['passed'] else 'FAILED'}"
              f" ({end_info['n_failed']} bin(s) differing by >= {END_MAX_PCT:g}% and >= {END_MIN_ABS:g}x)",
              f"  Largest bin difference    : {worst[5]:+.1f}% (positions {worst[0]}-{worst[1]})",
              f"  Step from position N to 1 : {end_info['step_f']}x folded, {end_info['step_o']}x rotated mapping",
              f"  Reported profile          : {profile_note}", ""]
    else:
        L += [f"End-window check (validation) : {profile_note}", ""]

    # ---------------------------- draft text ------------------------------
    rot_text = (", then rotated to begin at [start gene, e.g. cox1], and the rotated "
                "assembly was confirmed to be an exact cyclic shift of the original")
    if rotated == "yes":
        rot_clause = rot_text
    elif rotated == "no":
        rot_clause = ""
    else:                     # not recorded: leave the clause for the author
        rot_clause = "[" + rot_text + " -- keep only if rotated]"
    methods = (
        f"Read coverage was assessed with a circularity-aware mapping pipeline "
        f"(version {version}; {cite_in}). The assembly "
        f"was confirmed to contain no redundant terminal overlap from "
        f"circularization{rot_clause}. Reads were mapped with Bowtie2 "
        f"v{P('bowtie2', '[version]')} in local mode (Langmead & Salzberg 2012) to the assembly padded by "
        f"appending its first K bp to its 3' end (K = mean insert size + 4 SD, "
        f"minimum 300 bp, estimated with SAMtools v{P('samtools', '[version]')}; Danecek et al. "
        f"2021), so "
        f"that read pairs spanning the arbitrary origin of the linear sequence "
        f"could align. Only properly paired alignments were retained, with no "
        f"mapping-quality threshold, because reads within the duplicated window "
        f"receive low mapping quality by design. Depth from the appended region was "
        f"added back to the corresponding initial positions before coverage "
        f"statistics were calculated with NumPy (Harris et al. 2020) and pandas "
        f"(McKinney 2010), and figures were drawn with Matplotlib (Hunter 2007). The correction was checked against an unpadded mapping "
        f"across the genome interior, and against a second unpadded mapping, to the "
        f"assembly rotated by half its length, across the first and last K bp.")
    if end_info is not None and not end_info["passed"]:
        methods += (" Because the two differed by at least 10% and 3x in at least one "
                    "200-bp bin, depth in those end windows was taken from the "
                    "rotated mapping.")
    results = (f"The padding length (K) was {k} bp. Reads covered {100 * breadth / n:.1f}% "
               f"of the {n}-bp mitogenome, with a mean depth of {mean:.1f} +/- {sd:.1f}x "
               f"(median {median:.1f}x; range {dmin}-{dmax}x).")
    if interior is not None:
        pct, same = interior
        if abs(pct) < 0.0005 and same:
            results += (f" Across the genome interior (positions {k + 1}-{n - k}), the "
                        f"corrected and unpadded profiles were identical in mean and "
                        f"minimum depth, with minima at the same positions.")
        else:
            results += (f" Across the genome interior (positions {k + 1}-{n - k}), mean "
                        f"depth of the corrected and unpadded profiles differed by "
                        f"{abs(pct):.2f}%" + (", with minima at the same positions."
                                             if same else "."))
    if end_info is not None:
        w = end_info["w"]
        if end_info["passed"]:
            results += (f" Across the first and last {w} bp, the corrected profile "
                        f"agreed with the rotated mapping in every 200-bp bin (largest "
                        f"difference {abs(end_info['worst'][5]):.1f}%).")
        else:
            results += (f" Across the first and last {w} bp, the corrected profile "
                        f"departed from the rotated mapping in "
                        f"{end_info['n_failed']} of the 200-bp bins (largest difference "
                        f"{abs(end_info['worst'][5]):.1f}%), so depth at positions "
                        f"{n - w + 1}-{n} and 1-{w} was taken from the rotated mapping.")

    import textwrap
    L += ["DRAFT TEXT (check against the values above)",
          "-" * 72, "Methods:", *textwrap.wrap(methods, 72), "",
          "Results:", *textwrap.wrap(results, 72), ""]
    L += ["REFERENCES (for the draft text)",
          "-" * 72]
    for r in sorted(REFERENCES + [cite_ref]):
        L += textwrap.wrap(r, 72, subsequent_indent="    ", break_on_hyphens=False)
    L += ["", "Also cite the assembler you used, for example:"]
    for r in ASSEMBLER_REFERENCES:
        L += textwrap.wrap(r, 72, subsequent_indent="    ", break_on_hyphens=False)
    L += [""]
    return L


def write_profile_outputs(prefix, profile, n, k, header, meta, depth_path,
                          metadata_path, extra_inputs, profile_note):
    positions = np.arange(1, n + 1)

    folded_name = f"{prefix}_folded_depth.txt"
    folded_out_path = os.path.join(OUT_DIR, folded_name)
    pd.DataFrame({"Pos": positions, "Depth": profile}).to_csv(
        folded_out_path, sep="\t", index=False)

    mean = float(np.mean(profile))
    sd = float(np.std(profile))
    dmin = int(np.min(profile))
    dmax = int(np.max(profile))
    breadth = int(np.sum(profile > 0))

    stats_out_path = os.path.join(OUT_DIR, f"{prefix}_coverage_stats.txt")
    with open(stats_out_path, "w") as fh:
        fh.write(f"Reference: {header}\n")
        fh.write(f"Genome length (N): {n} bp\n")
        fh.write(f"Padding used (K): {k} bp\n")
        fh.write(f"Mean Depth: {mean:.4f}\n")
        fh.write(f"Std Dev: {sd:.4f}\n")
        fh.write(f"Min: {dmin}\n")
        fh.write(f"Max: {dmax}\n")
        fh.write(f"Breadth: {breadth} bp ({100 * breadth / n:.2f}% of genome)\n")
        fh.write(f"End-window check: {profile_note}\n")
        fh.write("\n".join(provenance_lines(
            "02_fold_and_plot.py", [depth_path, metadata_path] + extra_inputs,
            meta)) + "\n")

    print(f"Mean Depth: {mean:.2f}")
    print(f"Std Dev: {sd:.2f}")
    print(f"Min: {dmin}")
    print(f"Max: {dmax}")
    print(f"Breadth: {breadth} bp / {n} bp")
    print(f"End-window check: {profile_note}")
    print(f"Depth table -> {folded_out_path}")
    print(f"Stats -> {stats_out_path}")

    # =========================== LINEAR PLOT ================================
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.fill_between(positions, profile, color=LINEAR_COLOR, alpha=0.3)
    ax.plot(positions, profile, color=LINEAR_COLOR, lw=0.8)
    ax.axhline(y=mean, color=MEAN_LINE_COLOR, linestyle="--", lw=1.2,
               label=f"Mean: {mean:.1f}×")
    ax.set_title(PLOT_TITLE, fontsize=16)
    ax.set_xlabel("Genomic Position (bp)", fontsize=11)
    ax.set_ylabel("Read Depth", fontsize=11)
    ax.legend(fontsize=10)
    ax.set_xlim(positions.min(), positions.max())
    plt.tight_layout()
    linear_out_path = os.path.join(OUT_DIR, f"{prefix}_linear_coverage.png")
    plt.savefig(linear_out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Linear plot -> {linear_out_path}")

    # ========================== CIRCULAR PLOT ================================
    theta = 2 * math.pi * positions / n
    fig2 = plt.figure(figsize=(7, 7))
    ax2 = fig2.add_subplot(111, projection="polar")
    ax2.plot(theta, profile, color=LINEAR_COLOR, lw=0.8)
    ax2.fill(theta, profile, color=LINEAR_COLOR, alpha=0.3)
    mean_circle = np.full_like(theta, mean)
    ax2.plot(theta, mean_circle, color=MEAN_LINE_COLOR, linestyle="--", lw=1.2,
             label=f"Mean: {mean:.1f}×")
    ax2.set_theta_zero_location("N")
    ax2.set_theta_direction(-1)
    ax2.set_title(f"{PLOT_TITLE} (circular)", fontsize=12, pad=20)
    ax2.set_xlabel("")
    ax2.legend(loc="upper right", bbox_to_anchor=(1.15, 1.1), fontsize=9)
    plt.tight_layout()
    circular_out_path = os.path.join(OUT_DIR, f"{prefix}_circular_coverage.png")
    plt.savefig(circular_out_path, dpi=300, bbox_inches="tight")
    plt.close(fig2)
    print(f"Circular plot -> {circular_out_path}")
    return folded_name


def main():
    global IN_DIR, OUT_DIR, PLOT_TITLE
    ap = argparse.ArgumentParser(
        description="Fold padded depth, run the end-window and interior checks, "
                    "and write statistics and figures. With no arguments, reads "
                    "and writes next to this script and asks for a prefix.")
    ap.add_argument("--input-dir", help="folder holding the coverage_run/ files")
    ap.add_argument("--output-dir", help="folder for the outputs (created if needed)")
    ap.add_argument("--prefix", help="output prefix (otherwise asked)")
    ap.add_argument("--title", help="plot title (otherwise PLOT_TITLE below)")
    args = ap.parse_args()
    if args.input_dir:
        IN_DIR = os.path.abspath(args.input_dir)
    if args.output_dir:
        OUT_DIR = os.path.abspath(args.output_dir)
        os.makedirs(OUT_DIR, exist_ok=True)
    if args.title:
        PLOT_TITLE = args.title

    depth_path = os.path.join(IN_DIR, DEPTH_FILENAME)
    metadata_path = os.path.join(IN_DIR, METADATA_FILENAME)
    control_path = os.path.join(IN_DIR, CONTROL_FILENAME)
    origin_path = os.path.join(IN_DIR, ORIGIN_FILENAME)

    missing = [p for p in (depth_path, metadata_path) if not os.path.isfile(p)]
    if missing:
        print("ERROR: could not find the following required file(s) in")
        print(f"{IN_DIR}:")
        for p in missing:
            print(f"  - {os.path.basename(p)}")
        print("Copy padded_depth.txt and padding_metadata.txt into this folder")
        print("(from the coverage_run/ folder) and run this again.")
        return 1

    if args.prefix:
        out_prefix = args.prefix.strip()
    else:
        out_prefix = input(
            "Enter an output prefix for this sample (e.g. MySpecies_sample01): "
        ).strip()
    if not out_prefix:
        out_prefix = "coverage_output"
        print(f"No prefix entered -- using default: '{out_prefix}'")

    meta = read_metadata(metadata_path)
    n = int(meta["original_length_N"])
    k = int(meta["padding_length_K"])
    header = meta.get("original_header", "assembly")

    folded = fold_depth(read_depth_column(depth_path), n, k)

    # ======================= END-WINDOW CHECK (4.2+) =========================
    profile = folded
    extra_inputs = []
    end_lines = None
    end_info = None
    if not os.path.isfile(origin_path):
        profile_note = (f"not run ({ORIGIN_FILENAME} not found; produced by "
                        f"run_pipeline.sh version 4.2 or later)")
    elif "origin_check_shift_S" not in meta:
        profile_note = ("not run (padding_metadata.txt has no "
                        "origin_check_shift_S line; use the files from one run)")
    else:
        s = int(meta["origin_check_shift_S"])
        result = end_window_check(folded, read_depth_column(origin_path), n, k, s)
        if result is None:
            profile_note = "not run (see note printed on screen)"
        else:
            end_lines, passed, corrected, end_info = result
            extra_inputs = [origin_path]
            if passed:
                profile_note = "PASSED (folded profile reported)"
            else:
                profile_note = ("FAILED (end windows taken from the origin-check "
                                "mapping; see end-window check file)")
                profile = corrected
                keep = os.path.join(
                    OUT_DIR,
                    f"{out_prefix}_folded_depth_before_end_window_correction.txt")
                pd.DataFrame({"Pos": np.arange(1, n + 1), "Depth": folded}).to_csv(
                    keep, sep="\t", index=False)
                print(f"Folded profile kept for the record -> {keep}")

    folded_name = write_profile_outputs(out_prefix, profile, n, k, header, meta,
                                        depth_path, metadata_path, extra_inputs,
                                        profile_note)

    if end_lines is not None:
        end_out_path = os.path.join(OUT_DIR, f"{out_prefix}_end_window_check.txt")
        with open(end_out_path, "w") as fh:
            fh.write("\n".join(end_lines) + "\n")
            fh.write("\n".join(provenance_lines(
                "02_fold_and_plot.py", [depth_path, metadata_path, origin_path],
                meta)) + "\n")
        print("")
        print("\n".join(end_lines))
        print(f"End-window check -> {end_out_path}")
    else:
        print("")
        print(f"NOTE: end-window check {profile_note}.")

    # ===================== INTERIOR CHECK (VALIDATION) =======================
    # The end-window correction never touches the interior, so this compares
    # the reported profile and the folded profile alike.
    print("")
    control = None
    if not os.path.isfile(control_path):
        print(f"NOTE: {CONTROL_FILENAME} not found in {IN_DIR}, so the")
        print("interior check was skipped. To run it, copy that file there")
        print("(from the coverage_run/ folder) and run this script again.")
    else:
        control = read_depth_column(control_path)
        lines = compare_to_control(folded, control, n, k, folded_name)
        if lines is not None:
            comparison_out_path = os.path.join(OUT_DIR, f"{out_prefix}_control_comparison.txt")
            with open(comparison_out_path, "w") as fh:
                fh.write("\n".join(lines) + "\n")
                fh.write("\n".join(provenance_lines("02_fold_and_plot.py",
                                                     [depth_path, metadata_path, control_path],
                                                     meta)) + "\n")
            print("\n".join(lines))
            print(f"Interior check -> {comparison_out_path}")

    # ============================ REPORT ====================================
    prov_path = os.path.join(IN_DIR, PROVENANCE_FILENAME)
    report = build_report(out_prefix, profile, folded, control, n, k, meta,
                          read_provenance(prov_path), end_info, profile_note)
    inputs = [p for p in (depth_path, metadata_path, control_path, origin_path, prov_path)
              if os.path.isfile(p)]
    report_path = os.path.join(OUT_DIR, f"{out_prefix}_report.txt")
    with open(report_path, "w") as fh:
        fh.write("\n".join(report) + "\n")
        fh.write("\n".join(provenance_lines("02_fold_and_plot.py", inputs, meta)) + "\n")
    print("")
    print(f"Report (all values for the Methods and Results) -> {report_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
