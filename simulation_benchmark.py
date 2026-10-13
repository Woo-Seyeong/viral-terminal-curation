#!/usr/bin/env python3
"""
Simulation benchmark for the hairpin-terminal curation workflow (v3.1).

Ground truth is STRUCTURAL, not generative: the expected boundary at each
terminus is the end of the MAXIMALLY base-paired terminal stem. Bases at the
edge of the designed loop that happen to pair with their opposite number are
part of the stem by definition, so the workflow is scored against the maximal
stem, not against the loop as it was written. Scoring against the generative
loop would charge the workflow 1-2 bp of "error" for behaving correctly.

Usage
-----
    python simulation_benchmark.py /path/to/main.py [options]

`main.py` is the workflow code cell of hairpin_terminal_curator_v3.1.ipynb
saved as a .py file. Everything else is self-contained.

    --out DIR        output directory (default ./benchmark_out)
    --reps N         replicates per condition (default 3)
    --core-len N     internal genome length in bp (default 60000)
    --jobs N         parallel worker processes (default 4)
    --quick          small, fast grid for smoke-testing
    --conditions ... run only the named condition families

Outputs
-------
    results.csv      one row per simulated genome
    summary.txt      sensitivity, specificity and boundary error per condition
"""

import argparse
import csv
import io
import json
import contextlib
import multiprocessing as mp
import os
import random
import shutil
import sys
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


# --------------------------------------------------------------------------
# sequence construction
# --------------------------------------------------------------------------
def rc(s):
    return str(Seq(s).reverse_complement())


def rand_seq(n, gc, rng):
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ACGT", weights=w, k=n))


def at_rich_loop(n, rng, at_frac=0.9, edge_pairs=0):
    """An AT-rich terminal loop with a KNOWN amount of edge base pairing.

    Terminal loops of hairpin-ended dsDNA viruses are AT-rich, which makes
    chance complementarity between the two edges of a randomly generated loop
    common. That would leave the ground truth ambiguous, because such bases
    are part of the stem by definition. The loop is therefore built so that
    exactly `edge_pairs` positions at each edge pair with their opposite
    number and the next position does not: the expected stem boundary is then
    the designed boundary plus `edge_pairs`, exactly.
    """
    if n < 2 * edge_pairs + 3:
        raise ValueError("loop too short for the requested edge pairing")
    w = [at_frac / 2, (1 - at_frac) / 2, (1 - at_frac) / 2, at_frac / 2]
    loop = rng.choices("ACGT", weights=w, k=n)
    for k in range(edge_pairs):
        loop[n - 1 - k] = COMPLEMENT[loop[k]]
    k = edge_pairs
    while COMPLEMENT[loop[k]] == loop[n - 1 - k]:          # force the break
        loop[n - 1 - k] = rng.choice([b for b in "ACGT" if b != COMPLEMENT[loop[k]]])
    return "".join(loop)


def paired(a, b):
    return COMPLEMENT.get(a) == b


def maximal_stem(seq, loop_start, loop_len):
    """Extend perfect Watson-Crick pairing inward from a designed loop.

    Returns (stem_end, true_loop_len): stem_end is the first position that is
    not part of the maximally paired stem, i.e. the expected cut length.
    """
    k = 0
    while (
        k < loop_len // 2                      # the loop can shrink, not invert
        and paired(seq[loop_start + k], seq[loop_start + loop_len - 1 - k])
    ):
        k += 1
    # pairing extends k bases further in from each side, shrinking the loop
    return loop_start + k, loop_len - 2 * k


def mutate(seq, error_rate, rng):
    """Substitution noise, as a crude stand-in for residual basecall error."""
    if error_rate <= 0:
        return seq
    out = []
    for b in seq:
        if rng.random() < error_rate:
            out.append(rng.choice([x for x in "ACGT" if x != b]))
        else:
            out.append(b)
    return "".join(out)


def add_homopolymer_indel(seq, rng, n_events=1, max_len=2, min_len=1):
    """Shorten homopolymer runs, the ONT error mode that breaks the stem walk.

    This is the failure observed on the 3' terminus of ASFV HK_NT_202103,
    where a 9-T run opposed a 7-A run and the 20-bp stem window lost register.
    """
    import re

    runs = [m for m in re.finditer(r"([ACGT])\1{4,}", seq)]
    if not runs:
        return seq, 0
    out = seq
    for m in rng.sample(runs, min(n_events, len(runs))):
        d = rng.randint(min_len, max_len)
        s, e = m.start(), m.end()
        if e - s - d < 3:
            continue
        out = out[: e - d] + out[e:]
        return out, d  # one event keeps the coordinate bookkeeping simple
    return seq, 0


def insert_ns(seq, frac, rng):
    if frac <= 0:
        return seq
    n = int(len(seq) * frac)
    pos = rng.sample(range(len(seq)), n)
    out = list(seq)
    for p in pos:
        out[p] = "N"
    return "".join(out)


def build_genome(
    rng,
    core_len=60000,
    loop_len=35,
    loop_edge_pairs=0,
    gc=0.40,
    drc_left=0,
    drc_right=0,
    itr_len=0,
    error_rate=0.0,
    n_fraction=0.0,
    homopolymer_indel=False,
    homopolymer_indel_bp=0,
    internal_inverted_repeat=0,
):
    """Build a hairpin-ended genome and, optionally, its DRC-inflated version.

    Layout of the correctly represented ("true") genome:

        loopL + core + loopR

    where the sequence flanking each loop forms the terminal stem because the
    right end of the genome is, by construction, the reverse complement of the
    left end over the ITR.
    """
    body = rand_seq(core_len, gc, rng)

    if itr_len > 0:
        itr = rand_seq(itr_len, gc, rng)
        body = itr + body[itr_len:-itr_len] + rc(itr) if core_len > 2 * itr_len else body

    if internal_inverted_repeat > 0:
        # a decoy: an inverted repeat pair placed well inside the genome
        rep = rand_seq(internal_inverted_repeat, gc, rng)
        a = core_len // 4
        b = core_len // 2
        body = (
            body[:a]
            + rep
            + body[a + internal_inverted_repeat : b]
            + rc(rep)
            + body[b + internal_inverted_repeat :]
        )

    loopL = at_rich_loop(loop_len, rng, edge_pairs=loop_edge_pairs)
    loopR = at_rich_loop(loop_len, rng, edge_pairs=loop_edge_pairs)
    truth = loopL + body + loopR

    # A read that runs around the terminal hairpin continues into the
    # complementary strand, so the assembler emits the reverse complement of
    # the sequence INSIDE the loop. The loop itself stays at the junction and
    # is not duplicated.
    seq = truth
    if drc_left:
        seq = rc(body[:drc_left]) + seq
    if drc_right:
        seq = seq + rc(body[-drc_right:])

    # Expected (structural) boundaries, measured on the sequence as built:
    # the stem may extend into the designed loop wherever the edge bases pair.
    exp_left, true_loopL = maximal_stem(seq, drc_left, loop_len)
    rc_seq = rc(seq)
    exp_right, true_loopR = maximal_stem(rc_seq, drc_right, loop_len)

    # A homopolymer deletion inside the left DRC arm shortens everything that
    # precedes the left boundary, so the expected left cut moves by the number
    # of bases removed. The deletion is confined to the arm (never past the
    # boundary) so that the right boundary and the curated sequence are
    # unaffected and the comparison stays exact.
    if homopolymer_indel or homopolymer_indel_bp:
        region = min(400, max(0, exp_left - 50))
        d = homopolymer_indel_bp or 2
        head, deleted = add_homopolymer_indel(seq[:region], rng, max_len=d, min_len=d)
        seq = head + seq[region:]
        exp_left -= deleted
    seq = mutate(seq, error_rate, rng)
    seq = insert_ns(seq, n_fraction, rng)

    return {
        "seq": seq,
        "truth_len": len(truth),
        "expected_left_cut": exp_left,
        "expected_right_cut": exp_right,
        "designed_loop_len": loop_len,
        "structural_loop_len_L": true_loopL,
        "structural_loop_len_R": true_loopR,
    }


# --------------------------------------------------------------------------
# running the workflow
# --------------------------------------------------------------------------
_NS = {}


def load_workflow(path):
    ns = {}
    with open(path) as fh:
        code = fh.read()
    exec(code, ns)
    if "run_optimized_batch_pipeline" not in ns:
        raise SystemExit(
            f"{path} does not define run_optimized_batch_pipeline(); "
            "save the workflow code cell of the notebook as a .py file."
        )
    return ns


def run_one(args):
    workflow_path, tag, rec, params = args
    ns = _NS.get("ns")
    if ns is None:
        ns = load_workflow(workflow_path)
        _NS["ns"] = ns

    work = tempfile.mkdtemp(prefix="bench_")
    try:
        fa = os.path.join(work, "in.fa")
        with open(fa, "w") as fh:
            fh.write(f">{tag}\n{rec['seq']}\n")
        out = os.path.join(work, "out")
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(
            io.StringIO()
        ):
            ns["run_optimized_batch_pipeline"](
                fa,
                out,
                reference_fasta=None,
                min_itr_identity=0.95,
                orientation_min_ratio=1.1,
                itr_terminal_tolerance=500,
            )
        with open(os.path.join(out, "summary_report.csv")) as fh:
            row = next(csv.DictReader(fh))
        final = ""
        fp = os.path.join(out, "all_final_sequences.fasta")
        if os.path.exists(fp):
            seqs = list(SeqIO.parse(fp, "fasta"))
            if seqs:
                final = str(seqs[0].seq)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    i = lambda k: int(float(row.get(k) or 0))
    called_left = i("DRC_L_Cut")
    called_right = i("DRC_R_Cut")
    exp_l = rec["expected_left_cut"]
    exp_r = rec["expected_right_cut"]
    has_drc = (params["drc_left"] > 0) or (params["drc_right"] > 0)
    trimmed = (called_left + called_right) > 0

    return {
        **{k: params[k] for k in params},
        "rep_tag": tag,
        "input_len": len(rec["seq"]),
        "truth_len": rec["truth_len"],
        "final_len": i("Final_Length"),
        "trim_class": row.get("Trim_Class", ""),
        "expected_left_cut": exp_l,
        "expected_right_cut": exp_r,
        "called_left_cut": called_left,
        "called_right_cut": called_right,
        "phase2_left": i("L_STEP1_DRC_Cut_Bp"),
        "phase2_right": i("R_STEP1_DRC_Cut_Bp"),
        "phase3_left": i("L_STEP2_Refinement_Cut_Bp"),
        "phase3_right": i("R_STEP2_Refinement_Cut_Bp"),
        "left_error_bp": called_left - exp_l,
        "right_error_bp": called_right - exp_r,
        "designed_loop_len": rec["designed_loop_len"],
        "structural_loop_len_L": rec["structural_loop_len_L"],
        "structural_loop_len_R": rec["structural_loop_len_R"],
        "called_loop_len_L": i("L_Loop_Len"),
        "called_loop_len_R": i("R_Loop_Len"),
        "refine_method_L": row.get("L_Refine_Method", ""),
        "refine_method_R": row.get("R_Refine_Method", ""),
        "itr_status": row.get("ITR_Status", ""),
        "has_drc": int(has_drc),
        "trimmed": int(trimmed),
        "exact_sequence_recovered": int(
            final != "" and (final == rec["seq"][exp_l : len(rec["seq"]) - exp_r])
        ),
    }


# --------------------------------------------------------------------------
# condition grid
# --------------------------------------------------------------------------
def build_conditions(quick=False, wanted=None):
    """Each entry is (family, dict of build_genome kwargs)."""
    conds = []

    drc_lengths = [50, 100, 200, 300, 500, 1000, 2000, 5000, 10000, 20000, 50000]
    if quick:
        drc_lengths = [100, 500, 5000]
    for d in drc_lengths:
        conds.append(("drc_length", dict(drc_left=d, drc_right=d)))

    # one-sided artifacts: 4 of the public records were trimmed at the left
    # terminus only and 14 at the right only
    for d in ([500, 5000] if quick else [200, 1000, 5000, 20000]):
        conds.append(("one_sided_left", dict(drc_left=d, drc_right=0)))
        conds.append(("one_sided_right", dict(drc_left=0, drc_right=d)))

    for lp in ([20, 44] if quick else [11, 20, 25, 33, 44, 60]):
        conds.append(("loop_length", dict(drc_left=5000, drc_right=5000, loop_len=lp)))

    for gc in ([0.40] if quick else [0.30, 0.40, 0.50, 0.60]):
        conds.append(("gc_content", dict(drc_left=5000, drc_right=5000, gc=gc)))

    for er in ([0.0, 0.01] if quick else [0.0, 0.002, 0.005, 0.01, 0.02, 0.05]):
        conds.append(("error_rate", dict(drc_left=5000, drc_right=5000, error_rate=er)))

    for nf in ([0.0, 0.05] if quick else [0.0, 0.01, 0.05, 0.10, 0.25]):
        conds.append(("ambiguous_bases", dict(drc_left=5000, drc_right=5000, n_fraction=nf)))

    for itr in ([0, 2200] if quick else [0, 500, 2200, 5000]):
        conds.append(("itr_length", dict(drc_left=5000, drc_right=5000, itr_len=itr)))

    for ep in ([0, 3] if quick else [0, 1, 2, 3, 5]):
        conds.append(
            ("loop_edge_pairing", dict(drc_left=5000, drc_right=5000, loop_edge_pairs=ep))
        )

    # How large a length difference between the two stem arms the 20-bp stem
    # walk can absorb. The ASFV record in Section 3.8 differs by 2 bp at the
    # 3' terminus, so the ladder brackets that value.
    for d in ([2, 5] if quick else [1, 2, 3, 4, 5, 6, 8]):
        conds.append(
            ("homopolymer_indel", dict(drc_left=5000, drc_right=5000, homopolymer_indel_bp=d))
        )

    # --- negative controls: any trimming here is a false positive ---
    conds.append(("negative_clean", dict(drc_left=0, drc_right=0)))
    for rep in ([2000] if quick else [500, 2000, 10000]):
        conds.append(
            ("negative_internal_repeat", dict(drc_left=0, drc_right=0, internal_inverted_repeat=rep))
        )
    conds.append(("negative_high_n", dict(drc_left=0, drc_right=0, n_fraction=0.10)))
    conds.append(("negative_error", dict(drc_left=0, drc_right=0, error_rate=0.02)))

    if wanted:
        conds = [c for c in conds if c[0] in wanted]
    return conds


def summarize(rows, fh):
    from collections import defaultdict
    import statistics as st

    by = defaultdict(list)
    for r in rows:
        by[r["family"]].append(r)

    def p(*a):
        print(*a, file=fh)

    p("=" * 78)
    p("SUMMARY")
    p("=" * 78)
    pos = [r for r in rows if r["has_drc"]]
    neg = [r for r in rows if not r["has_drc"]]
    if pos:
        p(f"sensitivity (DRC present, any trimming): {sum(r['trimmed'] for r in pos)}/{len(pos)}")
    if neg:
        fp = sum(r["trimmed"] for r in neg)
        p(f"specificity (no DRC, untrimmed):         {len(neg) - fp}/{len(neg)}  false positives: {fp}")
    errs = [abs(r["left_error_bp"]) for r in pos] + [abs(r["right_error_bp"]) for r in pos]
    if errs:
        p(
            f"boundary error |bp| over positives: median {st.median(errs):.0f}  "
            f"mean {st.mean(errs):.1f}  max {max(errs)}  within 2 bp: "
            f"{sum(1 for e in errs if e <= 2)}/{len(errs)}"
        )
    exact = sum(r["exact_sequence_recovered"] for r in pos)
    if pos:
        p(f"exact sequence recovered (positives):    {exact}/{len(pos)}")
    p("")
    p(f"{'family':26s} {'n':>4s} {'detected':>9s} {'falsePos':>9s} {'medianErr':>10s} {'maxErr':>7s} {'exact':>9s}")
    for fam in sorted(by):
        rs = by[fam]
        pos = [r for r in rs if r["has_drc"]]
        neg = [r for r in rs if not r["has_drc"]]
        e = [abs(r["left_error_bp"]) for r in pos] + [abs(r["right_error_bp"]) for r in pos]
        med = f"{st.median(e):.0f}" if e else "-"
        mx = f"{max(e)}" if e else "-"
        det = f"{sum(r['trimmed'] for r in pos)}/{len(pos)}" if pos else "-"
        fp = f"{sum(r['trimmed'] for r in neg)}/{len(neg)}" if neg else "-"
        ex = f"{sum(r['exact_sequence_recovered'] for r in pos)}/{len(pos)}" if pos else "-"
        p(f"{fam:26s} {len(rs):4d} {det:>9s} {fp:>9s} {med:>10s} {mx:>7s} {ex:>9s}")
    hp = [r for r in rows if r["family"] == "homopolymer_indel"]
    if hp:
        p("")
        p("Homopolymer length difference between stem arms:")
        for d in sorted({r["homopolymer_indel_bp"] for r in hp}):
            g = [r for r in hp if r["homopolymer_indel_bp"] == d]
            e = [abs(r["left_error_bp"]) for r in g]
            p(f"  {d} bp: detected {sum(r['trimmed'] for r in g)}/{len(g)}  "
              f"left-boundary error median {st.median(e):.0f}, max {max(e)}")
    p("")
    p("Per-condition detail (DRC length family):")
    for r in sorted(
        (r for r in rows if r["family"] == "drc_length"), key=lambda r: r["drc_left"]
    ):
        p(
            f"  DRC={r['drc_left']:>6d} bp  class={r['trim_class']:<16s} "
            f"err L/R = {r['left_error_bp']:+d}/{r['right_error_bp']:+d}  "
            f"phase2 L/R = {r['phase2_left']}/{r['phase2_right']}"
        )


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("workflow", help="path to the workflow code saved as .py")
    ap.add_argument("--out", default="benchmark_out")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--core-len", type=int, default=60000)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20260101)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--conditions", nargs="*", default=None)
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    conds = build_conditions(quick=a.quick, wanted=a.conditions)

    jobs = []
    for fam, kw in conds:
        for rep in range(a.reps):
            rng = random.Random(f"{a.seed}:{fam}:{json.dumps(kw, sort_keys=True)}:{rep}")
            full = dict(
                core_len=a.core_len,
                loop_len=35,
                gc=0.40,
                drc_left=0,
                drc_right=0,
                itr_len=0,
                error_rate=0.0,
                n_fraction=0.0,
                homopolymer_indel=False,
                homopolymer_indel_bp=0,
                internal_inverted_repeat=0,
                loop_edge_pairs=0,
            )
            full.update(kw)
            rec = build_genome(rng, **full)
            params = dict(full)
            params["family"] = fam
            params["rep"] = rep
            tag = f"{fam}_{rep}_" + "_".join(f"{k}{v}" for k, v in sorted(kw.items()))
            jobs.append((a.workflow, tag[:90], rec, params))

    print(f"{len(conds)} conditions x {a.reps} replicates = {len(jobs)} genomes", file=sys.stderr)

    rows = []
    if a.jobs > 1:
        with mp.Pool(a.jobs) as pool:
            for n, r in enumerate(pool.imap_unordered(run_one, jobs), 1):
                rows.append(r)
                print(f"  {n}/{len(jobs)}", end="\r", file=sys.stderr)
    else:
        for n, j in enumerate(jobs, 1):
            rows.append(run_one(j))
            print(f"  {n}/{len(jobs)}", end="\r", file=sys.stderr)
    print(file=sys.stderr)

    cols = list(rows[0].keys())
    with open(os.path.join(a.out, "results.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(a.out, "summary.txt"), "w") as fh:
        summarize(rows, fh)
    summarize(rows, sys.stdout)
    print(f"\nwrote {a.out}/results.csv and {a.out}/summary.txt")


if __name__ == "__main__":
    main()
