#!/usr/bin/env python3
"""
All-vs-all TM-score within a single directory of structures, then cluster them.

Runs US-align in batch mode (parallelized by chunking the file list), records the
raw pairwise output to a file, parses the TM-scores, and clusters the structures
with a choice of graph/hierarchical methods.

Usage:
    python3 cluster_dir.py STRUCT_DIR                         # MCL, threshold 0.45
    python3 cluster_dir.py STRUCT_DIR --method mcl --inflation 1.6
    python3 cluster_dir.py STRUCT_DIR --method complete --threshold 0.45

The pairwise scores are saved to <dir>_pairs.tsv and REUSED automatically on later
runs (so re-clustering with a different --method/--inflation/--threshold is
instant). Pass --force to recompute, or --from-raw FILE to point at another
<dir>_pairs.tsv file.

Clustering methods (--method):
    mcl       (default) Markov Clustering on the TM-weighted graph (edges with
              TM > threshold, weighted by TM). Flow-based community detection;
              granularity controlled by --inflation (LOWER inflation -> fewer,
              LARGER clusters; higher -> more, smaller). Avoids single-linkage
              chaining while still merging genuinely related structures.
    single    connected components of the >threshold graph (single-linkage).
    complete  agglomerative complete-linkage cut at the threshold (conservative:
              every pair within a cluster exceeds the threshold).
    average   agglomerative average-linkage cut at the threshold.
    greedy    CD-HIT/qTMclust style representative assignment.

Representative per cluster (--rep):
    longest   (default) longest chain in the cluster (most complete model)
    central   member with the highest mean TM-score to the rest of the cluster
    attractor the MCL attractor (MCL only; falls back to 'central' otherwise)

TM-score normalization (for weights & threshold):
    --norm shorter (default) | longer | avg   ("shorter" = normalized by the
    shorter of the two chains; most sensitive.)

Outputs (current directory):
    <dir>_pairs.tsv      pairwise scores (structure1, structure2, TM_used, TM1, TM2,
                         L1, L2) -- the recorded result, reused automatically and
                         via --from-raw
    clusters.tsv         one row per structure: structure, cluster_id, representative
    cluster_reps.txt     one representative structure per cluster

Requires: US-align; and for --method mcl/complete/average: pip install numpy scipy
"""
import argparse, os, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor


def list_structures(d, suffix):
    names = sorted(f[:-len(suffix)] for f in os.listdir(d) if f.endswith(suffix))
    if not names:
        sys.exit(f"No '*{suffix}' files found in {d}")
    return names


def run_chunk(usalign, d, chunk_file, full_file, suffix, mol, infmt):
    cmd = [usalign, "-dir1", d + os.sep, chunk_file,
           "-dir2", d + os.sep, full_file,
           "-suffix", suffix, "-outfmt", "2"]
    if mol != "auto":
        cmd += ["-mol", mol]
    if infmt is not None:
        cmd += ["-infmt1", infmt, "-infmt2", infmt]
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        if not line or line.startswith("#") or line.startswith("PDBchain1"):
            continue
        rows.append(line)
    return rows


def base(cid):
    """'name.pdb:A' or 'dir/name.pdb:A' -> 'name' (structure id)."""
    left = cid.split(":")[0]
    fn = left.rsplit("/", 1)[-1]
    for ext in (".pdb", ".cif", ".ent"):
        if fn.endswith(ext):
            return fn[:-len(ext)]
    return fn


def pick_tm(tm1, tm2, l1, l2, norm):
    if norm == "shorter":
        return tm1 if l1 <= l2 else tm2
    if norm == "longer":
        return tm1 if l1 >= l2 else tm2
    return (tm1 + tm2) / 2.0


class UnionFind:
    def __init__(self, items):
        self.p = {x: x for x in items}
    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


# ---------------------------------------------------------------- clustering ---
def cluster_single(names, edges):
    uf = UnionFind(names)
    for a, b in edges:
        uf.union(a, b)
    comp = {}
    for n in names:
        comp.setdefault(uf.find(n), []).append(n)
    return list(comp.values())


def cluster_greedy(names, edges, lengths):
    adj = {n: set() for n in names}
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    order = sorted(names, key=lambda n: lengths.get(n, 0), reverse=True)
    assigned, reps = {}, []
    for n in order:
        placed = False
        for r in reps:
            if r in adj[n]:
                assigned[n] = r; placed = True; break
        if not placed:
            reps.append(n); assigned[n] = n
    clusters = {}
    for n, r in assigned.items():
        clusters.setdefault(r, []).append(n)
    return list(clusters.values())


def cluster_hierarchical(names, best, threshold, method):
    import numpy as np
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    n = len(names)
    if n == 1:
        return [[names[0]]]
    idx = {nm: i for i, nm in enumerate(names)}
    D = np.ones((n, n)); np.fill_diagonal(D, 0.0)
    for (a, b), v in best.items():
        i, j = idx[a], idx[b]
        d = 1.0 - v[0]
        D[i, j] = D[j, i] = d
    Z = linkage(squareform(D, checks=False), method=method)
    labels = fcluster(Z, t=1.0 - threshold, criterion="distance")
    clusters = {}
    for k, lab in enumerate(labels):
        clusters.setdefault(int(lab), []).append(names[k])
    return list(clusters.values())


def cluster_mcl(names, best, threshold, inflation,
                expansion=2, max_iter=300, prune=1e-4, tol=1e-6):
    import numpy as np
    n = len(names)
    if n == 1:
        return [[names[0]]]
    idx = {nm: i for i, nm in enumerate(names)}
    A = np.zeros((n, n), dtype=float)
    for (a, b), v in best.items():
        if v[0] > threshold:
            i, j = idx[a], idx[b]
            A[i, j] = A[j, i] = v[0]
    di = np.arange(n)
    A[di, di] = np.maximum(A.max(axis=1), 1e-6)      # self-loops stabilise MCL

    def normcols(M):
        s = M.sum(axis=0); s[s == 0] = 1.0
        return M / s

    M = normcols(A)
    for _ in range(max_iter):
        last = M
        M = M @ M                       # expansion (power = 2)
        M = np.power(M, inflation)      # inflation
        M = normcols(M)
        M[M < prune] = 0.0              # prune small values
        M = normcols(M)
        if np.abs(M - last).max() < tol:
            print("converged!")
            break
    attractors = [i for i in range(n) if M[i, i] > 1e-8] or list(range(n))
    clusters = {}
    for j in range(n):
        col = M[attractors, j]
        a = attractors[int(np.argmax(col))] if col.max() > 0 else j
        clusters.setdefault(a, []).append(names[j])
    members = list(clusters.values())
    attr = [names[a] for a in clusters.keys()]   # MCL attractor per cluster
    return members, attr


# ------------------------------------------------------ representative pick ---
def tm_between(best, a, b):
    if a == b:
        return 1.0
    v = best.get((a, b) if a < b else (b, a))
    return v[0] if v else 0.0


def choose_rep(members, rep_mode, lengths, best, attractor=None):
    """Pick a cluster representative.
        longest   : longest chain (most complete model)
        central   : member with the highest mean TM-score to the rest of the cluster
        attractor : the MCL attractor; falls back to 'central' if unavailable
    """
    if len(members) == 1:
        return members[0]
    if rep_mode == "attractor" and attractor is not None:
        return attractor
    if rep_mode == "longest":
        return max(members, key=lambda n: lengths.get(n, 0))
    # 'central' (also the fallback for 'attractor' when none was provided)
    def centrality(m):
        others = [x for x in members if x != m]
        return sum(tm_between(best, m, x) for x in others) / max(1, len(others))
    return max(members, key=centrality)


# --------------------------------------------------------------------- main ---
def parse_rows(rows, norm):
    best, lengths = {}, {}
    for line in rows:
        p = line.split("\t")
        if len(p) < 11:
            continue
        a, b = base(p[0]), base(p[1])
        if a == b:
            continue
        try:
            tm1, tm2 = float(p[2]), float(p[3])
            l1, l2 = int(p[8]), int(p[9])
        except ValueError:
            continue
        lengths[a] = l1; lengths[b] = l2
        tm = pick_tm(tm1, tm2, l1, l2, norm)
        key = (a, b) if a < b else (b, a)
        prev = best.get(key)
        if prev is None or tm > prev[0]:
            best[key] = (tm, tm1, tm2, l1, l2)
    return best, lengths


def parse_pairs_file(path, norm):
    """Read a saved <dir>_pairs.tsv (structure1, structure2, TM_used, TM1, TM2,
    L1, L2). TM_used is recomputed from TM1/TM2/L1/L2 with the current --norm."""
    best, lengths = {}, {}
    with open(path) as f:
        f.readline()  # header
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            a, b = p[0], p[1]
            try:
                tm1, tm2 = float(p[3]), float(p[4])
                l1, l2 = int(p[5]), int(p[6])
            except ValueError:
                continue
            lengths[a] = l1; lengths[b] = l2
            tm = pick_tm(tm1, tm2, l1, l2, norm)
            key = (a, b) if a < b else (b, a)
            prev = best.get(key)
            if prev is None or tm > prev[0]:
                best[key] = (tm, tm1, tm2, l1, l2)
    return best, lengths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("directory")
    ap.add_argument("--threshold", type=float, default=0.45)
    ap.add_argument("--method", default="mcl",
                    choices=["mcl", "single", "complete", "average", "greedy"])
    ap.add_argument("--linkage", default=None,
                    help="(deprecated alias for --method)")
    ap.add_argument("--inflation", type=float, default=2.0,
                    help="MCL inflation; lower -> fewer/larger clusters (default 2.0)")
    ap.add_argument("--rep", default="longest",
                    choices=["longest", "central", "attractor"],
                    help="cluster representative: longest chain (default), most "
                         "central by TM, or the MCL attractor (falls back to central)")
    ap.add_argument("--suffix", default=".pdb",
                    help="structure file extension, e.g. .pdb or .cif")
    ap.add_argument("--mol", default="RNA", choices=["RNA", "protein", "auto"])
    ap.add_argument("--infmt", default="auto", choices=["auto", "PDB", "mmCIF"])
    ap.add_argument("--norm", default="shorter", choices=["shorter", "longer", "avg"])
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--usalign", default="USalign")
    ap.add_argument("--from-raw", default=None,
                    help="explicit <dir>_pairs.tsv file to read scores from (skips US-align)")
    ap.add_argument("--force", action="store_true",
                    help="re-run US-align even if a saved <dir>_pairs.tsv exists")
    args = ap.parse_args()
    method = args.linkage or args.method

    d = args.directory.rstrip(os.sep)
    infmt = {"auto": None, "PDB": "0", "mmCIF": "3"}[args.infmt]
    names = list_structures(d, args.suffix)
    tag = os.path.basename(os.path.abspath(d))       # resolves "." to the real name
    pairs_file = f"{tag}_pairs.tsv"
    print(f"{len(names)} structures in {d}/  (mol={args.mol}, threshold={args.threshold}, "
          f"norm={args.norm}, method={method}"
          + (f", inflation={args.inflation}" if method == "mcl" else "")
          + f", rep={args.rep})")

    # ---- gather pairwise scores ----
    # Reuse the saved pairs file automatically: an explicit --from-raw, else the
    # default <dir>_pairs.tsv if it exists (unless --force). Otherwise run US-align.
    reuse = args.from_raw or (pairs_file if (os.path.exists(pairs_file) and not args.force)
                              else None)
    if reuse:
        best, lengths = parse_pairs_file(reuse, args.norm)
        print(f"reusing {len(best)} pairs from {reuse} (US-align skipped; "
              f"use --force to recompute)")
    else:
        tmpdir = tempfile.mkdtemp(prefix="cluster_")
        full_file = os.path.join(tmpdir, "all.lst")
        with open(full_file, "w") as f:
            f.write("\n".join(names) + "\n")
        njobs = max(1, min(args.jobs, len(names)))
        chunk_files = []
        for i in range(njobs):
            chunk = names[i::njobs]
            if not chunk:
                continue
            cf = os.path.join(tmpdir, f"chunk_{i}.lst")
            with open(cf, "w") as f:
                f.write("\n".join(chunk) + "\n")
            chunk_files.append(cf)
        rows = []
        with ThreadPoolExecutor(max_workers=njobs) as ex:
            futs = [ex.submit(run_chunk, args.usalign, d, cf, full_file,
                              args.suffix, args.mol, infmt) for cf in chunk_files]
            for fut in futs:
                rows.extend(fut.result())
        best, lengths = parse_rows(rows, args.norm)
        with open(pairs_file, "w") as f:         # record the result for reuse
            f.write("structure1\tstructure2\tTM_used\tTM1\tTM2\tL1\tL2\n")
            for (a, b), (tm, tm1, tm2, l1, l2) in sorted(best.items()):
                f.write(f"{a}\t{b}\t{tm:.4f}\t{tm1:.4f}\t{tm2:.4f}\t{l1}\t{l2}\n")
        print(f"ran US-align on {len(names)} structures ({args.jobs} jobs); "
              f"scores -> {pairs_file}")

    edges = [(a, b) for (a, b), v in best.items() if v[0] > args.threshold]

    # ---- cluster ----
    attr = None
    if method == "single":
        clusters = cluster_single(names, edges)
    elif method == "greedy":
        clusters = cluster_greedy(names, edges, lengths)
    elif method in ("complete", "average"):
        clusters = cluster_hierarchical(names, best, args.threshold, method)
    else:  # mcl
        clusters, attr = cluster_mcl(names, best, args.threshold, args.inflation)
    if attr is None:
        attr = [None] * len(clusters)

    # sort clusters by size (largest first), keeping attractors aligned
    order = sorted(range(len(clusters)), key=lambda i: len(clusters[i]), reverse=True)
    clusters = [clusters[i] for i in order]
    attr = [attr[i] for i in order]

    with open("clusters.tsv", "w") as f, open("cluster_reps.txt", "w") as fr:
        f.write("structure\tcluster_id\trepresentative\n")
        for cid, (members, at) in enumerate(zip(clusters, attr), 1):
            rep = choose_rep(members, args.rep, lengths, best, at)
            fr.write(rep + "\n")
            for m in sorted(members):
                f.write(f"{m}\t{cid}\t{rep}\n")

    sizes = sorted((len(c) for c in clusters), reverse=True)
    print(f"pairs: {len(best)} | edges > {args.threshold}: {len(edges)}")
    print(f"clusters: {len(clusters)} | largest: {sizes[0] if sizes else 0} | "
          f"singletons: {sum(1 for s in sizes if s == 1)}")
    print(f"wrote clusters.tsv, cluster_reps.txt"
          + ("" if reuse else f", {pairs_file}"))


if __name__ == "__main__":
    main()
