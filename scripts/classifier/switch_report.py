#!/usr/bin/env python3
"""M8d: per-class evidence for switching the decision from the rule cascade to the learned classifier.

    switch_report.py --eval <eval.npz> [--curated <eval_curated.npz>] --batch <results dir>
                     [--review review.tsv] [--per-pair 30]

1. Held-out quality per class (saved `fit.py evaluate` predictions, gated): precision and recall of the
   model at each confidence band against the cascade's, on library families (held out by order) and
   on curated consensi (held out by cluster).
2. What a switch would change in a `classify` batch run with the model (`model_*` columns):
   Unknowns the cascade leaves unresolved that the model calls, and cascade TE calls it relabels.
3. A review sample of disagreements with their evidence (`--review`).

The switch policy (`policy`): the cascade's non-TE calls stand, and so do its calls backed by a protein hit
in their own class; otherwise the model decides at ≥ HYBRID (≥ UNKNOWN where the cascade made no call) unless
it calls Retroposon or nonTE.
"""
import argparse
import collections
import csv
import glob
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fit import CLASSES, gate, target  # noqa: E402

BANDS = [("high", 0.99), ("medium", 0.9), ("low", 0.7)]
HYBRID = 0.7      # where the cascade made a TE call
UNKNOWN = 0.9     # where it made none
KEEP_CASCADE = ("Retroposon", "nonTE")  # model classes that never override the cascade
CASCADE_CALLS = {"known_te", "chimeric"}
REVIEW_COLUMNS = ["best_subject", "best_bitscore", "nucl_best_subject", "nucl_bitscore", "nucl_coverage",
                  "domain_architecture", "length", "tail", "tsd_length", "copy_loci"]


def policy(prob, casc, prot_bits):
    """Final class per family (-1 = no call); `prob` gated, `prot_bits` best protein bits per class."""
    pred, conf = prob.argmax(1), prob.max(1)
    keep = np.isin(pred, [CLASSES.index(c) for c in KEEP_CASCADE])
    backed = (casc >= 0) & (prot_bits[np.arange(len(casc)), np.maximum(casc, 0)] > 0)
    take = (conf >= np.where(casc >= 0, HYBRID, UNKNOWN)) & ~keep & ~backed & (casc != CLASSES.index("nonTE"))
    return np.where(take, pred, casc)


def held_out(path, title):
    d = np.load(path)
    y, casc = d["y"], d["cascade"]
    p = gate(d["prob"], casc)
    pred, conf = p.argmax(1), p.max(1)
    print(f"\n== held-out quality: {title} ({len(y):,} labelled families)")
    print(f"{'class':<11}{'n':>8}{'cascade P/R':>16}" + "".join(f"{f'model ≥{t} P/R':>18}" for _, t in BANDS))
    for k, c in enumerate(CLASSES):
        t = y == k
        if not t.any():
            continue
        cells = [f"{((casc == k) & t).sum() / max(1, (casc == k).sum()):.3f}/{((casc == k) & t).sum() / t.sum():.3f}"]
        for _, thr in BANDS:
            s = (pred == k) & (conf >= thr)
            cells.append(f"{(s & t).sum() / max(1, s.sum()):.3f}/{(s & t).sum() / t.sum():.3f}")
        print(f"{c:<11}{t.sum():>8,}" + "".join(f"{x:>{16 if i == 0 else 18}}" for i, x in enumerate(cells)))
    # hybrid: the model's call at >= HYBRID, else the cascade's (which may abstain)
    variants = {"cascade": casc, "hybrid": np.where(conf >= HYBRID, pred, casc)}
    if "prot_bits" in d.files:
        variants["policy"] = policy(p, casc, d["prot_bits"])
    print(f"hybrid: the model at ≥{HYBRID}, else the cascade; policy: the hybrid, but cascade non-TE calls, "
          f"protein-backed cascade calls and model {'/'.join(KEEP_CASCADE)} calls keep the cascade, and the model "
          f"needs ≥{UNKNOWN} where the cascade made no call")
    for name, calls in variants.items():
        called = calls >= 0
        print(f"  {name:<8} calls {called.mean():.1%}, agreement {(calls[called] == y[called]).mean():.3f}")
    print(f"{'class':<11}" + "".join(f"{f'{n} P/R':>18}" for n in variants))
    for k, c in enumerate(CLASSES):
        t = y == k
        if t.any():
            print(f"{c:<11}" + "".join(
                f"{f'{((v == k) & t).sum() / max(1, (v == k).sum()):.3f}/{((v == k) & t).sum() / t.sum():.3f}':>18}"
                for v in variants.values()))


def batch_rows(results):
    """classification.tsv rows with the species and, from features.tsv, `prot_bits` per class."""
    for path in sorted(glob.glob(os.path.join(results, "*", "*", "classification.tsv"))):
        species = path.split(os.sep)[-3]
        bits = {}
        for f in csv.DictReader(open(os.path.join(os.path.dirname(path), "features.tsv")), delimiter="\t"):
            bits[f["family_id"]] = {c: float(f.get(f"prot_bits_{c}") or 0) for c in CLASSES}
        for r in csv.DictReader(open(path), delimiter="\t"):
            r["species"] = species
            r["prot_bits"] = bits.get(r["family_id"], {})
            yield r


def takes_model(r, cc):
    """`policy` for one batch row: does the model's call replace the cascade's?"""
    mc = r["model_class"]
    return (float(r["model_probability"]) >= (HYBRID if cc else UNKNOWN) and mc not in KEEP_CASCADE and cc != "nonTE"
            and not (cc and r["prot_bits"].get(cc, 0) > 0))


def cascade_class(r):
    if r["disposition"] in ("satellite", "simple_repeat", "low_complexity"):
        return "nonTE"
    t = target(r["label"])
    return t if t and r["disposition"] in CASCADE_CALLS | {"host_or_artefact"} else None


def switch_effects(results, review_path, per_pair, seed):
    unknown = collections.Counter()      # (cascade disposition, model class, band)
    relabel = collections.Counter()      # (cascade class, model class, band) for cascade calls
    total_unknown = 0
    samples = collections.defaultdict(list)
    rng = random.Random(seed)
    seen = collections.Counter()
    switched = collections.Counter()     # (cascade class or "none", model class, band) under the policy
    families = 0
    for r in batch_rows(results):
        families += 1
        if not r.get("model_class"):
            continue
        band = r["model_confidence"]
        mc = r["model_class"]
        cc = cascade_class(r)
        if takes_model(r, cc) and mc != cc:
            switched[(cc or "none", mc, band)] += 1
        if r["original_label"].startswith("Unknown"):
            total_unknown += 1
            if cc is None:
                unknown[(r["disposition"], mc, band)] += 1
        if cc and mc != cc and band in ("high", "medium"):
            relabel[(cc, mc, band)] += 1
            key = (cc, mc)
            seen[key] += 1
            # reservoir sample per class pair
            if len(samples[key]) < per_pair:
                samples[key].append(r)
            elif rng.random() < per_pair / seen[key]:
                samples[key][rng.randrange(per_pair)] = r
    print(f"\n== Unknown families the cascade does not call ({sum(unknown.values()):,} of {total_unknown:,} Unknowns)")
    print("model class at each band (gated non-TE excluded):")
    print(f"{'class':<11}" + "".join(f"{b:>10}" for b, _ in BANDS) + f"{'none':>10}")
    for c in CLASSES:
        row = [sum(n for (d, m, b), n in unknown.items() if m == c and b == band) for band, _ in BANDS]
        row.append(sum(n for (d, m, b), n in unknown.items() if m == c and b == "none"))
        print(f"{c:<11}" + "".join(f"{v:>10,}" for v in row))
    by_disp = collections.Counter()
    for (d, m, b), n in unknown.items():
        if b in ("high", "medium"):
            by_disp[d] += n
    print("resolved at medium or high, by cascade disposition:", dict(by_disp.most_common()))
    print("\n== cascade calls the model relabels (medium or high):")
    pairs = collections.Counter()
    for (cc, mc, b), n in relabel.items():
        pairs[(cc, mc)] += n
    for (cc, mc), n in pairs.most_common(20):
        hi = relabel[(cc, mc, "high")]
        print(f"  {cc:<11}-> {mc:<11}{n:>8,}  (high {hi:,})")
    print(f"\n== the policy's switches ({sum(switched.values()):,} of {families:,} families; "
          "cascade class or none -> model class, by band):")
    pairs = collections.Counter()
    for (cc, mc, b), n in switched.items():
        pairs[(cc, mc)] += n
    for (cc, mc), n in pairs.most_common(20):
        print(f"  {cc:<11}-> {mc:<11}{n:>8,}  " + "  ".join(f"{b} {switched[(cc, mc, b)]:,}" for b, _ in BANDS))
    if review_path:
        with open(review_path, "w") as fh:
            cols = ["species", "family_id", "original_label", "disposition", "label", "confidence",
                    "model_label", "model_probability", "model_confidence"] + REVIEW_COLUMNS
            print("\t".join(cols), file=fh)
            for (cc, mc), rows in sorted(samples.items(), key=lambda kv: -len(kv[1])):
                for r in rows:
                    print("\t".join(r.get(c, "") for c in cols), file=fh)
        print(f"\nreview sample ({per_pair} per class pair): {review_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--eval", help="fit.py evaluate --save output (library families)")
    p.add_argument("--curated", help="the matching _curated.npz")
    p.add_argument("--batch", help="classify batch results dir (<species>/<acc>/classification.tsv)")
    p.add_argument("--review", help="write a disagreement review sample here")
    p.add_argument("--per-pair", type=int, default=30)
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args()
    if a.eval:
        held_out(a.eval, "library families, held out by order")
    if a.curated:
        held_out(a.curated, "curated consensi, held out by cluster")
    if a.batch:
        switch_effects(a.batch, a.review, a.per_pair, a.seed)


if __name__ == "__main__":
    main()
