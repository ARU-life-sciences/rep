#!/usr/bin/env python3
"""Fit and evaluate the M8 top-level classifier from `rep classify` feature tables (ROADMAP M8b).

A multinomial logistic regression over the top-level TE classes plus nonTE, fitted on
`features.tsv` rows whose original label is known. Features are transformed from the
`rep.features.v1` columns (log1p for scores, counts and lengths; a missing indicator for
optional columns), then standardised. Evaluation holds out whole clades (orders), never
families, and compares the model with the rule cascade at the same coverage.

    fit.py evaluate --libraries <batch>/results --taxonomy <tsv> [--refdb-holdout <dir>]
    fit.py fit      --libraries <batch>/results --taxonomy <tsv> [--refdb-holdout <dir>] --out classifier.json

`--libraries` holds <species>/<accession>/features.tsv (RepeatClassifier labels: noisy).
`--refdb-holdout` holds fold<k>_<chunk>/features.tsv from curated consensi classified against
a reference with their own cluster held out (clean labels). `--taxonomy` is a TSV with
species, order, class and phylum columns, used for the clade splits.
"""
import argparse
import collections
import csv
import glob
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize

SCHEMA = "rep.features.v1"
ID_COLUMNS = ["family_id", "original_label", "label", "disposition", "confidence"]
CLASSES = ["DNA", "RC", "LTR", "LINE", "SINE", "PLE", "Retroposon", "nonTE"]
NON_TE = {"Satellite", "Simple_repeat", "Low_complexity", "rRNA", "tRNA", "snRNA", "scRNA",
          "srpRNA", "RNA", "Other", "ARTEFACT"}
# Heavy-tailed columns get log1p before standardising.
LOG_PREFIXES = ("prot_bits_", "hmm_bits_", "nucl_bits_")
LOG_COLUMNS = {"prot_unclassified_bits", "host_domains", "host_best_score", "length", "ltr_length",
               "tir_length", "tandem_period", "tandem_copies", "longest_orf_aa", "n_orfs",
               "copy_loci", "copy_full_length", "genome_bp", "tsd_length"}
# Dispositions where the rule cascade commits to its label, and those that mean non-TE.
CASCADE_CALLS = {"known_te", "chimeric", "host_or_artefact"}
CASCADE_NON_TE = {"satellite", "simple_repeat", "low_complexity"}


def target(label):
    """Top-level class of a label, as `canonical_class`/`is_non_te` in the Rust code; None if unknown."""
    if not label or label.startswith("Unknown"):
        return None
    if label.startswith("LINE/Penelope"):
        return "PLE"
    top = label.split("/")[0]
    if top in NON_TE:
        return "nonTE"
    return top if top in CLASSES else None


def taxonomy(path):
    """species (lowercase, underscores) -> (order, class, phylum)."""
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        key = r["species"].lower().replace(" ", "_")
        out[key] = (r["order"] or "NA", r["class"] or "NA", r["phylum"] or "NA")
    return out


def lineage(species, tax):
    words = species.split("_")
    for key in ("_".join(words[:3]), "_".join(words[:2])):
        if key in tax:
            return tax[key]
    genus = [v for k, v in tax.items() if k.split("_")[0] == words[0]]
    return genus[0] if genus else ("NA", "NA", "NA")


def read_tables(paths, keep):
    """Numeric features (NaN for empty) and the ID columns of rows where keep(original_label)."""
    header, num, ids, where = None, [], [], []
    for p in paths:
        with open(p) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            if header is None:
                header = h
                oi = header.index("original_label")
                id_idx = [header.index(c) for c in ID_COLUMNS]
                num_idx = [i for i, c in enumerate(header) if c not in ID_COLUMNS]
            elif h != header:
                sys.exit(f"{p}: header differs from the first table")
            for line in fh:
                r = line.rstrip("\n").split("\t")
                if not keep(r[oi]):
                    continue
                num.append([float(r[i]) if r[i] != "" else np.nan for i in num_idx])
                ids.append([r[i] for i in id_idx])
                where.append(p)
    columns = [header[i] for i in num_idx]
    return columns, np.array(num, dtype=np.float64), ids, where


class Transform:
    """Column transforms, fitted on training rows: log1p, missing indicators, standardisation."""

    def __init__(self, columns, num):
        self.columns = columns
        self.log = np.array([c.startswith(LOG_PREFIXES) or c in LOG_COLUMNS for c in columns])
        self.optional = [c for c, m in zip(columns, np.isnan(num).any(axis=0)) if m]
        self.opt_idx = [columns.index(c) for c in self.optional]
        x = self.raw(num)
        self.mean = x.mean(axis=0)
        self.std = x.std(axis=0)
        self.std[self.std == 0] = 1.0

    def names(self):
        return self.columns + [f"{c}_missing" for c in self.optional]

    def raw(self, num):
        missing = np.isnan(num[:, self.opt_idx]).astype(np.float64)
        x = np.nan_to_num(num, nan=0.0)
        x[:, self.log] = np.log1p(np.maximum(x[:, self.log], 0))
        return np.hstack([x, missing])

    def __call__(self, num):
        return (self.raw(num) - self.mean) / self.std

    def to_json(self):
        return {"columns": self.columns, "log1p": [c for c, lg in zip(self.columns, self.log) if lg],
                "missing_indicators": self.optional, "features": self.names(),
                "mean": self.mean.tolist(), "std": self.std.tolist()}


def softmax(z):
    z = z - z.max(axis=1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)


def fit_lr(x, y, weights, l2):
    """Weighted multinomial logistic regression (L-BFGS); returns (K, d+1) with the bias last."""
    n, d = x.shape
    k = len(CLASSES)
    xb = np.hstack([x, np.ones((n, 1))])
    onehot = np.eye(k)[y]
    wsum = weights.sum()

    def loss(flat):
        w = flat.reshape(k, d + 1)
        p = softmax(xb @ w.T)
        nll = -(weights * np.log(p[np.arange(n), y] + 1e-12)).sum() / wsum
        reg = 0.5 * l2 * (w[:, :d] ** 2).sum()
        grad = ((p - onehot) * weights[:, None]).T @ xb / wsum
        grad[:, :d] += l2 * w[:, :d]
        return nll + reg, grad.ravel()

    res = minimize(loss, np.zeros(k * (d + 1)), jac=True, method="L-BFGS-B",
                   options={"maxiter": 500})
    return res.x.reshape(k, d + 1)


def predict(w, x):
    return softmax(np.hstack([x, np.ones((len(x), 1))]) @ w.T)


def class_weights(y, power):
    """Inverse class frequency to `power` (0 = none, 1 = balanced), normalised to mean 1."""
    counts = np.bincount(y, minlength=len(CLASSES)).astype(float)
    per = np.where(counts > 0, (counts.sum() / np.maximum(counts, 1)) ** power, 0)
    w = per[y]
    return w / w.mean()


def load(args):
    """Labelled rows with their source, species and clade."""
    tax = taxonomy(args.taxonomy) if args.taxonomy else {}
    paths, meta = [], {}
    for p in sorted(glob.glob(os.path.join(args.libraries, "*", "*", "features.tsv"))) if args.libraries else []:
        species = p.split(os.sep)[-3]
        paths.append(p)
        meta[p] = ("library", species, lineage(species, tax))
    for p in sorted(glob.glob(os.path.join(args.refdb_holdout, "*", "features.tsv"))) if args.refdb_holdout else []:
        fold = p.split(os.sep)[-2].split("_")[0]
        paths.append(p)
        meta[p] = ("refdb", fold, ("refdb_" + fold, "refdb", "refdb"))
    columns, num, ids, where = read_tables(paths, lambda label: target(label) is not None)
    y = np.array([CLASSES.index(target(r[1])) for r in ids])
    info = [meta[p] for p in where]
    print(f"labelled rows: {len(ids):,} from {len(paths)} tables", file=sys.stderr)
    return columns, num, ids, y, info


def cascade_call(ids):
    """The cascade's top-level call per family, or -1 where it abstains."""
    li, di = ID_COLUMNS.index("label"), ID_COLUMNS.index("disposition")
    nonte = CLASSES.index("nonTE")

    def one(r):
        if r[di] in CASCADE_NON_TE:
            return nonte
        t = target(r[li])
        return CLASSES.index(t) if t and r[di] in CASCADE_CALLS else -1
    return np.array([one(r) for r in ids])


def clade_folds(info, k):
    """Assign orders to k folds greedily by size, so each fold holds out whole orders."""
    sizes = collections.Counter(order for _, _, (order, _, _) in info)
    load_ = [0] * k
    fold_of = {}
    for order, n in sizes.most_common():
        f = min(range(k), key=lambda i: load_[i])
        fold_of[order] = f
        load_[f] += n
    return np.array([fold_of[order] for _, _, (order, _, _) in info]), fold_of


def report(y, prob, casc, title):
    n = len(y)
    pred = prob.argmax(axis=1)
    conf = prob.max(axis=1)
    called = casc >= 0
    cov_c = called.mean()
    acc_c = (casc[called] == y[called]).mean() if called.any() else float("nan")
    # model at the cascade's coverage, and at the cascade's accuracy
    order = np.argsort(-conf)
    correct = (pred[order] == y[order]).cumsum() / np.arange(1, n + 1)
    m = max(1, int(round(cov_c * n)))
    acc_m = correct[m - 1]
    ok = np.where(correct >= acc_c)[0]
    cov_m = (ok.max() + 1) / n if len(ok) else 0.0
    print(f"\n== {title}: {n:,} labelled families")
    print(f"cascade: calls {cov_c:.1%}, agreement {acc_c:.3f}")
    print(f"model  : at the same coverage agreement {acc_m:.3f}; at the cascade's agreement coverage {cov_m:.1%};"
          f" all families {(pred == y).mean():.3f}")
    print(f"{'class':<11}{'n':>9}{'casc_rec':>9}{'casc_prec':>10}{'model_rec':>10}{'model_prec':>11}   (model at cascade coverage)")
    thr = conf[order[m - 1]]
    sel = conf >= thr
    for k, c in enumerate(CLASSES):
        t = y == k
        if not t.any():
            continue
        cr = (called & t & (casc == y)).sum() / t.sum()
        cp = ((casc == k) & t).sum() / max(1, (casc == k).sum())
        mr = (sel & t & (pred == y)).sum() / t.sum()
        mp = (sel & (pred == k) & t).sum() / max(1, (sel & (pred == k)).sum())
        print(f"{c:<11}{t.sum():>9,}{cr:>9.3f}{cp:>10.3f}{mr:>10.3f}{mp:>11.3f}")


def evaluate(args):
    columns, num, ids, y, info = load(args)
    casc = cascade_call(ids)
    src = np.array([s for s, _, _ in info])
    lib = src == "library"
    folds, fold_of = clade_folds([i for i, s in zip(info, lib) if s], args.folds)
    fold = np.full(len(y), -1)
    fold[lib] = folds
    prob = np.zeros((len(y), len(CLASSES)))
    for f in range(args.folds):
        test = fold == f
        train = ~test  # refdb rows (clean labels) are always in training here
        held = sorted(o for o, k in fold_of.items() if k == f)
        print(f"fold {f}: train {train.sum():,}, test {test.sum():,}; orders held out: {', '.join(held[:12])}"
              f"{' …' if len(held) > 12 else ''}", file=sys.stderr)
        xi = np.where(train)[0][:: args.subsample]
        tr = Transform(columns, num[xi])
        w = class_weights(y[xi], args.balance) * np.where(src[xi] == "refdb", args.refdb_weight, 1.0)
        coef = fit_lr(tr(num[xi]), y[xi], w, args.l2)
        ti = np.where(test)[0]
        prob[ti] = predict(coef, tr(num[ti]))
    report(y[lib], prob[lib], casc[lib], "library families, held out by order")
    if args.save:
        np.savez_compressed(args.save, y=y[lib], prob=prob[lib], cascade=casc[lib], classes=CLASSES,
                            family=[r[0] for r, l in zip(ids, lib) if l],
                            species=[sp for (s, sp, _), l in zip(info, lib) if l],
                            order=[o for (s, _, (o, _, _)), l in zip(info, lib) if l])
    if (~lib).any():
        refdb_eval(args, columns, num, y, info, casc)


def refdb_eval(args, columns, num, y, info, casc):
    """Curated consensi: train on libraries plus the other refdb folds, test on each refdb fold."""
    src = np.array([s for s, _, _ in info])
    rfold = np.array([sp if s == "refdb" else "" for s, sp, _ in info])
    prob = np.zeros((len(y), len(CLASSES)))
    for f in sorted(set(rfold) - {""}):
        test = rfold == f
        xi = np.where(~test)[0][:: args.subsample]
        tr = Transform(columns, num[xi])
        w = class_weights(y[xi], args.balance) * np.where(src[xi] == "refdb", args.refdb_weight, 1.0)
        coef = fit_lr(tr(num[xi]), y[xi], w, args.l2)
        ti = np.where(test)[0]
        prob[ti] = predict(coef, tr(num[ti]))
    r = src == "refdb"
    report(y[r], prob[r], casc[r], "curated consensi, held out by cluster fold")


def fit(args):
    columns, num, ids, y, info = load(args)
    src = np.array([s for s, _, _ in info])
    xi = np.arange(len(y))[:: args.subsample]
    tr = Transform(columns, num[xi])
    w = class_weights(y[xi], args.balance) * np.where(src[xi] == "refdb", args.refdb_weight, 1.0)
    coef = fit_lr(tr(num[xi]), y[xi], w, args.l2)
    out = {"schema": "rep.classifier.v1", "features_schema": SCHEMA, "model": "multinomial_logistic",
           "classes": CLASSES, "transform": tr.to_json(),
           "weights": coef[:, :-1].tolist(), "bias": coef[:, -1].tolist(),
           "training": {"rows": int(len(xi)), "per_class": {c: int((y[xi] == k).sum()) for k, c in enumerate(CLASSES)},
                        "l2": args.l2, "balance": args.balance, "refdb_weight": args.refdb_weight}}
    json.dump(out, open(args.out, "w"), indent=1)
    names = tr.names()
    for k, c in enumerate(CLASSES):
        top = np.argsort(-coef[k, :-1])[:8]
        print(f"{c:<11}" + ", ".join(f"{names[i]} {coef[k, i]:+.2f}" for i in top))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("command", choices=["evaluate", "fit"])
    p.add_argument("--libraries")
    p.add_argument("--refdb-holdout")
    p.add_argument("--taxonomy")
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--l2", type=float, default=1e-4)
    p.add_argument("--balance", type=float, default=0.5, help="class weight power (0 none, 1 balanced)")
    p.add_argument("--refdb-weight", type=float, default=1.0, help="weight of curated rows relative to library rows")
    p.add_argument("--subsample", type=int, default=1, help="use every n-th training row")
    p.add_argument("--out", default="classifier.json")
    p.add_argument("--save", help="evaluate: write held-out predictions to this .npz")
    args = p.parse_args()
    {"evaluate": evaluate, "fit": fit}[args.command](args)


if __name__ == "__main__":
    main()
