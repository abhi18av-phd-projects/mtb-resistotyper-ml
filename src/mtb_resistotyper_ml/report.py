"""Human-readable rendering of a prediction.

The structured result is the machine artefact; this is the laboratory-facing
view of the same object. It exists to make two things impossible to miss: which
tier the call came from, and which of the reasons are genuine drug mechanisms
rather than co-travellers.
"""

from __future__ import annotations


def render(result: dict) -> str:
    lines: list[str] = []
    add = lines.append
    conf = result.get("confidence", {})
    add(f"sample     {result['sample_id']}")
    add(f"drug       {result['drug']}")
    add(f"call       {result['prediction']}   P(R) = {result['p_resistant']}")

    tier, auc = conf.get("tier"), conf.get("auc")
    if tier and auc is not None:
        add(f"tier       {tier}  (AUC {auc:.3f}, nested lineage-leave-one-out)")
    else:
        add("tier       UNDECLARED — this bundle states no operating range")

    reasons = result.get("reasons", [])
    add("")
    add(f"reasons    {len(reasons)} mutation(s) carried")
    for r in reasons:
        target = "on-target" if r["on_target"] else "OFF-target"
        add(f"  {r['mutation']:<20} coef {r['coefficient']:+.3f}   "
            f"logit {r['logit_contribution']:+.3f}   {target:<10} {r['role']}")
    if not reasons:
        add("  none — the call rests on the intercept alone")

    add("")
    add("note       " + conf.get("caveat", ""))
    return "\n".join(lines)
