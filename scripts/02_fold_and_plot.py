#!/usr/bin/env python3
"""
02_fold_and_plot.py

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
     padding acts. Phase 7 of run_pipeline.sh therefore maps the reads to a
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
This version is written to be run directly from the Windows Python IDLE
shell (Run > Run Module, or F5) -- no command-line arguments needed.

Before running, put these files in the SAME FOLDER as this script (all
from the coverage_run/ folder written by run_pipeline.sh):
    padded_depth.txt            required
    padding_metadata.txt        required
    prelim_filtered_depth.txt   optional, needed for the interior check
    origin_check_depth.txt      optional, needed for the end-window check
When you run it, the program will ask you to type an output prefix (e.g.
MySpecies_sample01) -- this is used to name all the output files.

It only needs pandas, matplotlib, and numpy -- no bioinformatics tools --
so it also runs fine inside WSL if you'd rather use that instead.

Outputs (written next to this script, named using your output prefix)
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

03_compare_to_control.py produces the same interior comparison as a
standalone step, for depth profiles folded with an earlier version of this
script.

Provenance: the statistics and check files end with a provenance block
recording the pipeline version, date, Python and library versions, and MD5
checksums of the input files (and, via padding_metadata.txt, of the
assembly). The values above the block are unchanged by it.
"""

import hashlib
import os
import math
import platform
from datetime import datetime

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

__version__ = "4.2.2"   # coverage pipeline version

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

# Input files are expected in the same folder as this script.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEPTH_FILENAME = "padded_depth.txt"
METADATA_FILENAME = "padding_metadata.txt"
CONTROL_FILENAME = "prelim_filtered_depth.txt"
ORIGIN_FILENAME = "origin_check_depth.txt"


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
    return lines, passed, corrected


def write_profile_outputs(prefix, profile, n, k, header, meta, depth_path,
                          metadata_path, extra_inputs, profile_note):
    positions = np.arange(1, n + 1)

    folded_name = f"{prefix}_folded_depth.txt"
    folded_out_path = os.path.join(SCRIPT_DIR, folded_name)
    pd.DataFrame({"Pos": positions, "Depth": profile}).to_csv(
        folded_out_path, sep="\t", index=False)

    mean = float(np.mean(profile))
    sd = float(np.std(profile))
    dmin = int(np.min(profile))
    dmax = int(np.max(profile))
    breadth = int(np.sum(profile > 0))

    stats_out_path = os.path.join(SCRIPT_DIR, f"{prefix}_coverage_stats.txt")
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
    linear_out_path = os.path.join(SCRIPT_DIR, f"{prefix}_linear_coverage.png")
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
    circular_out_path = os.path.join(SCRIPT_DIR, f"{prefix}_circular_coverage.png")
    plt.savefig(circular_out_path, dpi=300, bbox_inches="tight")
    plt.close(fig2)
    print(f"Circular plot -> {circular_out_path}")
    return folded_name


def main():
    depth_path = os.path.join(SCRIPT_DIR, DEPTH_FILENAME)
    metadata_path = os.path.join(SCRIPT_DIR, METADATA_FILENAME)
    control_path = os.path.join(SCRIPT_DIR, CONTROL_FILENAME)
    origin_path = os.path.join(SCRIPT_DIR, ORIGIN_FILENAME)

    missing = [p for p in (depth_path, metadata_path) if not os.path.isfile(p)]
    if missing:
        print("ERROR: could not find the following required file(s) in this")
        print(f"script's folder ({SCRIPT_DIR}):")
        for p in missing:
            print(f"  - {os.path.basename(p)}")
        print("Copy padded_depth.txt and padding_metadata.txt into this folder")
        print("(from the coverage_run/ folder) and run this again.")
        return

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
            end_lines, passed, corrected = result
            extra_inputs = [origin_path]
            if passed:
                profile_note = "PASSED (folded profile reported)"
            else:
                profile_note = ("FAILED (end windows taken from the origin-check "
                                "mapping; see end-window check file)")
                profile = corrected
                keep = os.path.join(
                    SCRIPT_DIR,
                    f"{out_prefix}_folded_depth_before_end_window_correction.txt")
                pd.DataFrame({"Pos": np.arange(1, n + 1), "Depth": folded}).to_csv(
                    keep, sep="\t", index=False)
                print(f"Folded profile kept for the record -> {keep}")

    folded_name = write_profile_outputs(out_prefix, profile, n, k, header, meta,
                                        depth_path, metadata_path, extra_inputs,
                                        profile_note)

    if end_lines is not None:
        end_out_path = os.path.join(SCRIPT_DIR, f"{out_prefix}_end_window_check.txt")
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
    if not os.path.isfile(control_path):
        print(f"NOTE: {CONTROL_FILENAME} not found next to this script, so the")
        print("interior check was skipped. To run it, copy that file here")
        print("(from the coverage_run/ folder) and run this script again.")
        return

    control = read_depth_column(control_path)
    lines = compare_to_control(folded, control, n, k, folded_name)
    if lines is None:
        return
    comparison_out_path = os.path.join(SCRIPT_DIR, f"{out_prefix}_control_comparison.txt")
    with open(comparison_out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
        fh.write("\n".join(provenance_lines("02_fold_and_plot.py",
                                             [depth_path, metadata_path, control_path],
                                             meta)) + "\n")
    print("\n".join(lines))
    print(f"Interior check -> {comparison_out_path}")


if __name__ == "__main__":
    main()
