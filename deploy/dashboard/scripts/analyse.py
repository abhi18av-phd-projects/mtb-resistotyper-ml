#!/usr/bin/env python3
"""Turn the scored cohort into the artefact the dashboard reads.

Adds the two things the raw scores do not say on their own: a threshold-free AUC
on THIS cohort, and the ROC sweep that shows where the operating point should
have been. Both exist because the fixed 0.5 cut turned an informative ranking
into a near-zero sensitivity, and a dashboard that reported only the calls would
have shown a failure where there is a miscalibration.
"""
import json, sys
from pathlib import Path

R = Path(sys.argv[1])
rows = json.loads((R / "predictions.json").read_text())
summary = json.loads((R / "summary.json").read_text())

by = {}
for r in rows:
    by.setdefault(r["drug"], []).append(r)

def auc(pairs):
    pos = [p for p, y in pairs if y]; neg = [p for p, y in pairs if not y]
    if not pos or not neg: return None
    return round(sum((a > b) + 0.5 * (a == b) for a in pos for b in neg) / (len(pos) * len(neg)), 4)

for d, rs in by.items():
    pairs = [(r["p_resistant"], r["phenotype"] == "R") for r in rs]
    n_r = sum(1 for _, y in pairs if y); n_s = len(pairs) - n_r
    summary[d]["auc_on_unclassified"] = auc(pairs)
    summary[d]["prevalence"] = round(n_r / len(pairs), 4)
    roc = []
    for t in sorted({p for p, _ in pairs}):
        tp = sum(1 for p, y in pairs if p >= t and y)
        fp = sum(1 for p, y in pairs if p >= t and not y)
        roc.append({"threshold": t,
                    "sensitivity": round(tp / n_r, 4) if n_r else None,
                    "specificity": round((n_s - fp) / n_s, 4) if n_s else None,
                    "tp": tp, "fp": fp})
    summary[d]["roc"] = roc
    # The operating point a deployer would actually pick: the most sensitive cut
    # that still holds specificity at 0.90. Reported alongside the inherited 0.5
    # so the cost of the inherited one is legible.
    ok = [p for p in roc if p["specificity"] is not None and p["specificity"] >= 0.90]
    best = max(ok, key=lambda p: p["sensitivity"]) if ok else None
    summary[d]["recalibrated"] = best

(R / "summary.json").write_text(json.dumps(summary, indent=2))
tot = {k: sum(summary[d][k] for d in summary) for k in ("tp", "n_resistant")}
rec = sum((summary[d]["recalibrated"] or {}).get("tp", 0) for d in summary)
print(f"at the inherited 0.5 cut : {tot['tp']}/{tot['n_resistant']} resistant recovered")
print(f"at a recalibrated cut    : {rec}/{tot['n_resistant']} recovered (specificity held at 0.90)")
