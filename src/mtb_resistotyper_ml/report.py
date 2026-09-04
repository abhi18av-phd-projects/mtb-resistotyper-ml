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
    layer = result.get("layer", 2)
    add(f"sample     {result['sample_id']}")
    add(f"drug       {result['drug']}")
    add(f"source     {result.get('source', 'model')}  (layer {layer})")

    p_r = result.get("p_resistant")
    # A catalogue grade has no probability. Printing "P(R) = None" invites the
    # reader to treat the absence as a missing number rather than as a different
    # kind of answer.
    add(f"call       {result['prediction']}"
        + (f"   P(R) = {p_r}" if p_r is not None else "   graded, not modelled"))

    tier, auc = conf.get("tier"), conf.get("auc")
    if layer == 1:
        cat = result.get("provenance", {}).get("catalogue") or {}
        add(f"catalogue  {cat.get('catalogue', 'WHO')} {cat.get('version', '')}"
            f"  ({cat.get('grammar', 'GARC1')})")
    elif tier and auc is not None:
        add(f"tier       {tier}  (AUC {auc:.3f}, nested lineage-leave-one-out)")
    else:
        add("tier       UNDECLARED — this bundle states no operating range")

    rel = result.get("reliability")
    if rel:
        bar = "#" * rel["band"] + "." * (rel["of"] - rel["band"])
        capped = ("  (capped: intercept only)" if rel.get("capped_by_evidence")
                  else "  (capped by tier)" if rel["capped_by_tier"] else "")
        # A catalogue grade has no probability and therefore no margin. Its band
        # is 5 because a curated entry is not a close call, not because two
        # classes were far apart.
        margin = ("" if rel.get("margin") is None
                  else f"  margin {rel['margin']:.3f}")
        add(f"reliability {rel['band']}/{rel['of']}  [{bar}]{margin}{capped}")

    reasons = result.get("reasons", [])
    add("")
    add(f"reasons    {len(reasons)} mutation(s) carried")
    for r in reasons:
        target = "on-target" if r["on_target"] else "OFF-target"
        if r.get("logit_contribution") is None:
            # A catalogue-graded mutation has a call, not a coefficient. Printing
            # a zero here would suggest the catalogue weighed it and found it
            # negligible, which inverts what a graded entry means.
            add(f"  {r['mutation']:<20} WHO call {r.get('call', '?'):<3}"
                f"          {target:<10} {r['role']}")
        else:
            add(f"  {r['mutation']:<20} coef {r['coefficient']:+.3f}   "
                f"logit {r['logit_contribution']:+.3f}   {target:<10} {r['role']}")
    if not reasons:
        add("  none — the call rests on the intercept alone"
            if result.get("layer") != 1 else "  none graded by the catalogue")

    chain = result.get("causal_chain")
    if chain:
        from mtb_resistotyper_ml.explain import render_chain
        add("")
        add("reasoning chain")
        add(render_chain(chain))

    add("")
    add("note       " + conf.get("caveat", ""))
    return "\n".join(lines)
