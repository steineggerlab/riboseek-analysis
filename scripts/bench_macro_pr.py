#!/usr/bin/env python3
"""
bench_macro_pr.py - Macro-averaged precision/recall benchmark (Option 2).

TRUE per-query macro average (not the old normalized-ratio precision):

    For each query q, using whichever hits are "included" at the current
    point on the curve:
        TP_q = # included hits in the SAME family   (self-hit allowed -> TP)
        FP_q = # included hits in a DIFFERENT family (not clan-ignored)
        precision_q = TP_q / (TP_q + FP_q)
        recall_q    = TP_q / nF(q)          # nF(q) = size of q's family

    macro precision = mean over queries of precision_q
    macro recall    = mean over queries of recall_q

nF(q) (the number of possible true positives per query) is precomputed from
the lookup file.

Two ways to define "included" hits, selected with --mode:

  --mode evalue   (DEFAULT, recommended)
      Sweep a single GLOBAL score threshold t. Every query is judged at the
      same stringency, which is fairer across families of different sizes.
      macro precision is averaged over queries that have at least one
      prediction at t (N_PRED); recall is averaged over all eligible
      queries (N_QUERIES). At the loosest threshold both converge to the
      final macro numbers.

      The direction depends on --score:
        --score evalue (default): column 3 is an e-value, LOWER is better,
            a hit is included if e-value <= t.
        --score bits: column 3 is a bit score, HIGHER is better, a hit is
            included if bit-score >= t.

  --mode depth
      x-axis is per-query rank depth k = 1,2,3... (top-k hits per query),
      carry-forward once a query runs out of hits. Uses only per-query
      ordering (no global sort). Kept for comparison.

Rules (both modes):
  * Self hits allowed, count as TP.
  * FAMILY-level classification only.
  * If the lookup has a clan column (CLxxxx) and two ids are in different
    families but the same clan, the hit is IGNORED (not TP, not FP).
  * Optional --query-list restricts scoring to listed queries.

Usage:
    python3 bench_macro_pr.py --lookup LOOKUP.tsv --m8 HITS.m8 --out PREFIX \
        [--mode evalue|depth] [--query-list Q.list] [--step N]

Outputs PREFIX.tsv and PREFIX.png.
"""
import argparse
import sys
from collections import defaultdict, OrderedDict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_lookup(path):
    id2fam, id2clan, famCnt = {}, {}, defaultdict(int)
    has_clan = False
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] == "":
                continue
            sid, fam = parts[0], parts[1]
            clan = parts[2] if len(parts) >= 3 else ""
            id2fam[sid] = fam
            id2clan[sid] = clan
            famCnt[fam] += 1
            if clan.startswith("CL"):
                has_clan = True
    return id2fam, id2clan, famCnt, has_clan


def read_hits(args, id2fam, id2clan, famCnt, has_clan, qlist):
    """Return per-query list of (evalue, is_tp) for valid, non-ignored hits.

    Also returns counters. Ignored (same-clan/diff-family) hits are dropped.
    """
    perq = OrderedDict()          # q -> list of (evalue, is_tp)
    n_skipped_no_fam = 0
    n_ignored_clan = 0
    with open(args.m8) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            q, t, ev = p[0], p[1], p[2]
            if qlist is not None and q not in qlist:
                continue
            if q not in id2fam or t not in id2fam:
                n_skipped_no_fam += 1
                continue
            fam_q, fam_t = id2fam[q], id2fam[t]
            if fam_q == fam_t:
                is_tp = True
            elif has_clan:
                cq, ct = id2clan.get(q, ""), id2clan.get(t, "")
                if cq != "" and ct != "" and cq == ct:
                    n_ignored_clan += 1
                    continue
                is_tp = False
            else:
                is_tp = False
            try:
                evf = float(ev)
            except ValueError:
                continue
            perq.setdefault(q, []).append((evf, is_tp))
    return perq, n_skipped_no_fam, n_ignored_clan


def curve_depth(perq, id2fam, famCnt, higher_is_better=False):
    """Mean precision@k / recall@k over queries, carry-forward. Returns
    (xvals, macroP, macroR, nq)."""
    queries = list(perq.keys())
    nq = len(queries)
    # order each query's hits best-first (e-value ascending, or bits descending)
    key = (lambda x: -x[0]) if higher_is_better else (lambda x: x[0])
    prec_seq, rec_seq = {}, {}
    for q in queries:
        hits = sorted(perq[q], key=key)
        nF = max(famCnt.get(id2fam[q], 1), 1)
        tp = fp = 0
        ps, rs = [], []
        for ev, is_tp in hits:
            if is_tp:
                tp += 1
            else:
                fp += 1
            ps.append(tp / (tp + fp))
            rs.append(tp / nF)
        prec_seq[q], rec_seq[q] = ps, rs
    K = max(len(prec_seq[q]) for q in queries)
    sumP, sumR = np.zeros(K + 1), np.zeros(K + 1)
    for q in queries:
        ps, rs = prec_seq[q], rec_seq[q]
        m = len(ps)
        for k in range(1, m + 1):
            sumP[k] += ps[k - 1]
            sumR[k] += rs[k - 1]
        if m < K:
            sumP[m + 1:K + 1] += ps[-1]
            sumR[m + 1:K + 1] += rs[-1]
    xvals = np.arange(1, K + 1)
    return xvals, sumP[1:] / nq, sumR[1:] / nq, nq, None


def curve_threshold(perq, id2fam, famCnt, higher_is_better=False):
    """Sweep a global score threshold. At each distinct score emit a macro
    precision/recall point, processing hits best-first (e-value ascending, or
    bit-score descending). Returns (thresholds, macroP, macroR, nq, npred)."""
    queries = list(perq.keys())
    nq = len(queries)
    nF = {q: max(famCnt.get(id2fam[q], 1), 1) for q in queries}

    # flatten all hits, sort best-first
    all_hits = []
    for q in queries:
        for sc, is_tp in perq[q]:
            all_hits.append((sc, q, is_tp))
    all_hits.sort(key=(lambda x: -x[0]) if higher_is_better else (lambda x: x[0]))

    tp = defaultdict(int)
    fp = defaultdict(int)
    sum_prec = 0.0     # sum over queries-with-predictions of precision_q
    sum_rec = 0.0      # sum over all queries of recall_q
    n_pred = 0

    thr, Pc, Rc, NP = [], [], [], []
    i, N = 0, len(all_hits)
    while i < N:
        ev = all_hits[i][0]
        # process all hits at this exact e-value (a cutoff t includes e <= t)
        while i < N and all_hits[i][0] == ev:
            _, q, is_tp = all_hits[i]
            old_tot = tp[q] + fp[q]
            old_prec = (tp[q] / old_tot) if old_tot > 0 else 0.0
            old_rec = tp[q] / nF[q]
            if is_tp:
                tp[q] += 1
            else:
                fp[q] += 1
            new_tot = tp[q] + fp[q]
            new_prec = tp[q] / new_tot
            new_rec = tp[q] / nF[q]
            if old_tot == 0:
                n_pred += 1          # query gains its first prediction
            sum_prec += new_prec - old_prec
            sum_rec += new_rec - old_rec
            i += 1
        thr.append(ev)
        Pc.append(sum_prec / n_pred if n_pred > 0 else 0.0)
        Rc.append(sum_rec / nq)
        NP.append(n_pred)
    return np.array(thr), np.array(Pc), np.array(Rc), nq, np.array(NP)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lookup", required=True)
    ap.add_argument("--m8", required=True)
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--mode", choices=["evalue", "depth"], default="evalue",
                    help="'evalue' = global score-threshold sweep (default); "
                         "'depth' = per-query top-k")
    ap.add_argument("--score", choices=["evalue", "bits"], default="evalue",
                    help="type of score in column 3: 'evalue' (lower is "
                         "better, default) or 'bits' (higher is better)")
    ap.add_argument("--query-list", default=None)
    ap.add_argument("--step", type=int, default=1,
                    help="write every Nth curve point to the tsv (default 1)")
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    id2fam, id2clan, famCnt, has_clan = load_lookup(args.lookup)
    qlist = None
    if args.query_list:
        qlist = set(x.strip() for x in open(args.query_list) if x.strip())

    perq, n_skip, n_ign = read_hits(args, id2fam, id2clan, famCnt, has_clan, qlist)
    if not perq:
        sys.exit("No eligible queries with valid hits found - check inputs.")

    title = args.title or args.out
    higher_is_better = (args.score == "bits")
    score_label = "bit-score" if higher_is_better else "e-value"

    if args.mode == "evalue":
        xv, macroP, macroR, nq, npred = curve_threshold(
            perq, id2fam, famCnt, higher_is_better)
        xlabel_curve = f"{score_label} threshold"
        col = "BITS" if higher_is_better else "EVALUE"
        header = f"{col}\tMACRO_PREC\tMACRO_RECALL\tN_QUERIES\tN_PRED"
    else:
        xv, macroP, macroR, nq, npred = curve_depth(
            perq, id2fam, famCnt, higher_is_better)
        xlabel_curve = "per-query depth k (top-k hits)"
        header = "DEPTH\tMACRO_PREC\tMACRO_RECALL\tN_QUERIES"

    K = len(xv)
    tsv_path = args.out + ".tsv"
    with open(tsv_path, "w") as out:
        out.write(header + "\n")
        for j in range(K):
            if (j + 1) % args.step == 0 or j == K - 1:
                if args.mode == "evalue":
                    sval = (f"{xv[j]:.6g}" if higher_is_better else f"{xv[j]:.3e}")
                    out.write(f"{sval}\t{macroP[j]:.6f}\t{macroR[j]:.6f}\t{nq}\t{npred[j]}\n")
                else:
                    out.write(f"{int(xv[j])}\t{macroP[j]:.6f}\t{macroR[j]:.6f}\t{nq}\n")

    # plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].plot(macroR, macroP, marker="o", markersize=3, linewidth=1)
    axes[0].set_xlabel("Macro recall (mean per-query)")
    axes[0].set_ylabel("Macro precision (mean per-query)")
    axes[0].set_title(f"Precision vs Recall - {title}")
    axes[0].grid(True, alpha=0.3)

    if args.mode == "evalue" and not higher_is_better:
        # x = e-value threshold on log scale; guard zeros
        xpos = np.array(xv, dtype=float)
        tiny = xpos[xpos > 0].min() if np.any(xpos > 0) else 1e-300
        xpos = np.where(xpos <= 0, tiny, xpos)
        axes[1].plot(xpos, macroP, label="Macro precision", linewidth=1)
        axes[1].plot(xpos, macroR, label="Macro recall", linewidth=1)
        axes[1].set_xscale("log")
        axes[1].set_xlabel(xlabel_curve + " (log scale)")
    elif args.mode == "evalue" and higher_is_better:
        # bit-score threshold, linear; invert so strict(high) is on the left
        axes[1].plot(xv, macroP, label="Macro precision", linewidth=1)
        axes[1].plot(xv, macroR, label="Macro recall", linewidth=1)
        axes[1].invert_xaxis()
        axes[1].set_xlabel(xlabel_curve + " (strict -> loose)")
    else:
        axes[1].plot(xv, macroP, label="Macro precision", linewidth=1)
        axes[1].plot(xv, macroR, label="Macro recall", linewidth=1)
        axes[1].set_xlabel(xlabel_curve)
    axes[1].set_ylabel("Value")
    axes[1].set_title(f"Precision / Recall vs threshold - {title}")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    png_path = args.out + ".png"
    plt.savefig(png_path, dpi=150)

    print(f"[{title}] mode={args.mode}", file=sys.stderr)
    print(f"  clan column detected: {has_clan}", file=sys.stderr)
    print(f"  eligible queries: {nq}", file=sys.stderr)
    print(f"  curve points: {K}", file=sys.stderr)
    print(f"  hits skipped (id without family): {n_skip}", file=sys.stderr)
    print(f"  hits ignored (same clan, diff family): {n_ign}", file=sys.stderr)
    print(f"  FINAL (loosest) macro precision: {macroP[-1]:.6f}", file=sys.stderr)
    print(f"  FINAL (loosest) macro recall:    {macroR[-1]:.6f}", file=sys.stderr)
    print(f"  wrote {tsv_path} and {png_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
