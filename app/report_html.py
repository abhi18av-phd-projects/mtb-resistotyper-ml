"""One self-contained HTML report, in place of a folder of files to correlate.

The zip carried six artefacts and left the reader to join them: a TSV of calls,
a JSON of contributions, a text rendering, the chains, the inputs. Each is
correct and none answers the question a laboratory actually opens the results
to ask, which is "what is this isolate resistant to, and why should I believe
it". This assembles that answer, in the shape MultiQC established for
bioinformatics reporting: one file, offline, no assets to lose.

Design decisions worth stating because they are easy to get wrong:

* Resistance is a STATUS, not a category or a magnitude, so it uses the reserved
  status palette and always carries the letter alongside the colour. A reader
  with deuteranopia, a greyscale printout, or forced-colors reads the same grid.

* A contribution plot is DIVERGING, not sequential: the question is which side of
  zero a mutation pushed, and by how much. Blue and red with a neutral midpoint,
  never a one-hue ramp, which would make "pushes toward susceptible" look merely
  small rather than opposite.

* Catalogue and model are two identities, so they take categorical slots one and
  three, validated as a pair. They are also labelled in text everywhere, because
  the whole point of the report is that the reader can tell which layer answered.

* No CDN. The report travels in a zip and is opened later, offline, possibly on a
  machine with no route to the internet; a chart library fetched at open time
  would render a blank page exactly when it matters.
"""

from __future__ import annotations

import datetime
import html
import json
from typing import Any

# ── palette (validated: adjacent CVD ΔE 23.1, normal-vision 24.0, light) ──────
LIGHT = {
    "surface": "#fcfcfb", "plane": "#f9f9f7", "ink": "#0b0b0b", "ink2": "#52514e",
    "muted": "#898781", "grid": "#e1e0d9", "axis": "#c3c2b7",
    "border": "rgba(11,11,11,0.10)",
    "catalogue": "#2a78d6", "model": "#1baf7a",
    "good": "#0ca30c", "critical": "#d03b3b", "warning": "#fab219",
    "div_pos": "#e34948", "div_neg": "#2a78d6", "div_mid": "#f0efec",
}
DARK = {
    "surface": "#1a1a19", "plane": "#0d0d0d", "ink": "#ffffff", "ink2": "#c3c2b7",
    "muted": "#898781", "grid": "#2c2c2a", "axis": "#383835",
    "border": "rgba(255,255,255,0.10)",
    "catalogue": "#3987e5", "model": "#199e70",
    "good": "#0ca30c", "critical": "#d03b3b", "warning": "#fab219",
    "div_pos": "#e66767", "div_neg": "#3987e5", "div_mid": "#383835",
}


def _e(x: Any) -> str:
    return html.escape(str(x), quote=True)


def _vars(d: dict) -> str:
    return "".join(f"--{k}:{v};" for k, v in d.items())


def _waterfall(reasons: list[dict], width: int = 560) -> str:
    """Signed logit contributions as a diverging bar chart about zero.

    Sorted by magnitude rather than by name: the reader's question is which
    mutation moved the call, and rank answers it in one glance where alphabetical
    order answers nothing. Bars are labelled directly, so the diverging hues are
    reinforcement rather than the only channel.
    """
    rows = [r for r in reasons if r.get("logit_contribution") is not None]
    if not rows:
        return '<p class="empty">No modelled mutation was carried; the call rests on the intercept alone.</p>'
    rows = sorted(rows, key=lambda r: -abs(r["logit_contribution"]))
    top = max(abs(r["logit_contribution"]) for r in rows) or 1.0
    bar_h, gap = 20, 9
    # Size the label gutter to the longest name rather than fixing it, so a plot
    # of short gene names is not two-thirds empty.
    pad_l = min(210, max(96, int(max(len(r["mutation"]) for r in rows) * 6.6) + 16))
    plot = width - pad_l - 66
    mid = pad_l + plot / 2
    h = len(rows) * (bar_h + gap) + 26

    out = [f'<svg class="wf" viewBox="0 0 {width} {h}" role="img" '
           f'aria-label="Signed log-odds contribution of each carried mutation">']
    out.append(f'<line x1="{mid}" y1="4" x2="{mid}" y2="{h - 22}" '
               f'stroke="var(--axis)" stroke-width="1"/>')
    for i, r in enumerate(rows):
        y = i * (bar_h + gap) + 4
        v = r["logit_contribution"]
        w = max(2.0, abs(v) / top * (plot / 2 - 6))
        x = mid if v >= 0 else mid - w
        colour = "var(--div_pos)" if v >= 0 else "var(--div_neg)"
        on = r.get("on_target")
        out.append(
            f'<text x="{pad_l - 10}" y="{y + bar_h - 6}" text-anchor="end" '
            f'class="wfl {"on" if on else "off"}">{_e(r["mutation"])}</text>')
        # 4px rounded data-end anchored to the zero baseline, square at the axis.
        # An off-target contributor is drawn lighter as well as labelled: the
        # distinction between mechanism and co-occurrence is the report's main
        # claim, so it gets two channels, not one.
        dim = "" if on else ' fill-opacity="0.55"'
        coef = r.get("coefficient")
        tip = (f'{r["mutation"]} — {"on-target" if on else "off-target"}\n'
               f'coefficient {coef:+.3f}\nlog-odds {v:+.3f}\n{r.get("role", "")}')
        out.append(
            f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{bar_h}" rx="4" '
            f'fill="{colour}"{dim}><title>{_e(tip)}</title></rect>')
        out.append(f'<text x="{(x + w + 7) if v >= 0 else (x - 7):.1f}" y="{y + bar_h - 6}" '
                   f'text-anchor="{"start" if v >= 0 else "end"}" class="wfv">{v:+.2f}</text>')
    out.append(f'<text x="{mid - 8}" y="{h - 6}" text-anchor="end" class="wfa">'
               f'← toward susceptible</text>')
    out.append(f'<text x="{mid + 8}" y="{h - 6}" class="wfa">toward resistant →</text>')
    out.append("</svg>")
    return "".join(out)


def _meter(rel: dict) -> str:
    b, of = rel.get("band", 0), rel.get("of", 5)
    cap = ("capped: intercept only" if rel.get("capped_by_evidence")
           else "capped by tier" if rel.get("capped_by_tier") else "")
    seg = "".join(
        f'<span class="seg{" on" if i < b else ""}"></span>' for i in range(of))
    note = f' <span class="cap">{_e(cap)}</span>' if cap else ""
    return (f'<span class="meter" title="{_e(rel.get("basis", ""))}">{seg}'
            f'<span class="mv">{b}/{of}</span>{note}</span>')


def _chain(chain: list[dict]) -> str:
    if not chain:
        return ""
    out = ['<ol class="chain">']
    for s in chain:
        out.append(f'<li><span class="stage">{_e(s["stage"])}</span>'
                   f'<span class="kind">{_e(s["kind"])}</span>'
                   f'<p>{_e(s["claim"])}</p>')
        d = s.get("detail")
        if s["stage"] == "contribution" and isinstance(d, dict):
            for c in d.get("causal_reading_supported", []):
                out.append(f'<p class="mech"><b>{_e(c["mutation"])}</b> '
                           f'{c["logit"]:+.3f} — {_e(c["why"])}</p>')
            for c in d.get("association_only", []):
                out.append(f'<p class="assoc"><b>{_e(c["mutation"])}</b> '
                           f'{c["logit"]:+.3f} — {_e(c["why"])}</p>')
        elif isinstance(d, list) and d:
            out.append(f'<p class="det">{_e(", ".join(str(x) for x in d))}</p>')
        elif isinstance(d, dict):
            bits = [f"{k}: {v}" for k, v in d.items() if v not in (None, "", [], {})]
            if bits:
                out.append(f'<p class="det">{_e(" · ".join(bits))}</p>')
        out.append(f'<p class="basis">basis: {_e(s["basis"])}</p></li>')
    out.append("</ol>")
    return "".join(out)


def _antibiogram(by_sample: dict[str, list[dict]], drugs: list[str]) -> str:
    """Samples × drugs, the thing a laboratory reads first.

    A grid, not a chart: seven or more classes that all carry meaning belong in a
    table. Resistance is a status, so it gets the reserved status colours AND the
    letter in every cell. The layer is carried by a hairline underline plus the
    legend, never by hue alone, because a reader who cannot separate the two
    layers cannot tell a curated grade from a model estimate.
    """
    head = "".join(f"<th>{_e(d)}</th>" for d in drugs)
    rows = []
    for sid, rs in by_sample.items():
        by_drug = {r["drug"]: r for r in rs}
        cells = []
        for d in drugs:
            r = by_drug.get(d)
            if r is None:
                cells.append('<td class="ab none" title="not scored">·</td>')
                continue
            call = r["prediction"]
            layer = r.get("layer", 2)
            p = r.get("p_resistant")
            tip = (f'{sid} · {d}\n{"resistant" if call == "R" else "susceptible"}\n'
                   f'{"WHO catalogue (layer 1)" if layer == 1 else "model (layer 2)"}'
                   + (f'\nP(resistant) = {p}' if p is not None else ''))
            cells.append(
                f'<td class="ab {"r" if call == "R" else "s"} l{layer}">'
                f'<span title="{_e(tip)}">{call}</span></td>')
        rows.append(f'<tr><th class="sid">{_e(sid)}</th>{"".join(cells)}</tr>')
    return (f'<table class="antibiogram"><thead><tr><th></th>{head}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table>')


def _tiles(rows: list[dict]) -> str:
    """A KPI row, not a chart. Four headline numbers do not want a bar chart."""
    n_r = sum(1 for r in rows if r["prediction"] == "R")
    n_cat = sum(1 for r in rows if r.get("layer") == 1)
    n_mod = len(rows) - n_cat
    n_samples = len({r["sample_id"] for r in rows})
    tiles = [
        ("Isolates", n_samples, "", ""),
        ("Resistant calls", n_r, f"of {len(rows)}", "crit" if n_r else ""),
        ("Graded by catalogue", n_cat, "layer 1", "cat"),
        ("Answered by model", n_mod, "layer 2", "mod"),
    ]
    return '<div class="tiles">' + "".join(
        f'<div class="tile {cls}"><span class="lab">{_e(lab)}</span>'
        f'<span class="fig">{_e(v)}</span>'
        f'<span class="sub">{_e(sub)}</span></div>'
        for lab, v, sub, cls in tiles) + "</div>"


def _drug_card(r: dict) -> str:
    layer = r.get("layer", 2)
    conf = r.get("confidence", {})
    call = r["prediction"]
    badge = ("WHO catalogue" if layer == 1 else "model")
    head = (f'<summary><span class="dr">{_e(r["drug"])}</span>'
            f'<span class="call {"r" if call == "R" else "s"}">{call}</span>'
            f'<span class="lyr l{layer}">{_e(badge)} · L{layer}</span>')
    if r.get("p_resistant") is not None:
        head += f'<span class="pr">P(R) {r["p_resistant"]}</span>'
    if conf.get("auc") is not None:
        head += f'<span class="pr">AUC {conf["auc"]}</span>'
    head += _meter(r.get("reliability", {})) + "</summary>"

    body = []
    if layer == 1:
        graded = r.get("reasons", [])
        if graded:
            body.append('<table class="graded"><thead><tr><th>Mutation</th>'
                        '<th>WHO call</th></tr></thead><tbody>')
            body += [f'<tr><td>{_e(g["mutation"])}</td>'
                     f'<td class="{"r" if g.get("call") == "R" else "s"}">'
                     f'{_e(g.get("call", "?"))}</td></tr>' for g in graded]
            body.append("</tbody></table>")
        if r.get("ungraded_variants"):
            body.append(f'<p class="det">{len(r["ungraded_variants"])} further variant(s) '
                        f'were not gradable for this drug and did not affect the call.</p>')
    else:
        body.append(_waterfall(r.get("reasons", [])))
        if r.get("causal_chain"):
            body.append('<details class="chainwrap"><summary class="mini">'
                        'Reasoning chain — how this call was reached</summary>'
                        + _chain(r["causal_chain"]) + "</details>")
    body.append(f'<p class="caveat">{_e(conf.get("caveat", ""))}</p>')
    return f'<details class="card">{head}{"".join(body)}</details>'


CSS = """
*{box-sizing:border-box}
html{color-scheme:light}
body{margin:0;background:var(--plane);color:var(--ink);
     font:14.5px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}
.wrap{max-width:960px;margin:0 auto;padding:34px 22px 80px}
h1{font-size:24px;margin:0 0 3px;letter-spacing:-.01em}
h2{font-size:11.5px;text-transform:uppercase;letter-spacing:.08em;
   color:var(--muted);margin:34px 0 12px;font-weight:700}
.sub{color:var(--ink2);margin:0 0 22px;font-size:13.5px}
.panel{background:var(--surface);border:1px solid var(--border);border-radius:11px;
       padding:18px;margin-bottom:14px}
.det,.empty{color:var(--ink2);font-size:13px;margin:8px 0}
.caveat{color:var(--ink2);font-size:12.5px;margin:12px 0 0;padding-top:10px;
        border-top:1px solid var(--grid)}

/* KPI row — headline numbers, not a chart */
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:11px}
.tile{background:var(--surface);border:1px solid var(--border);border-radius:11px;padding:14px 16px}
.tile .lab{display:block;font-size:11px;text-transform:uppercase;letter-spacing:.07em;color:var(--muted)}
.tile .fig{display:block;font-size:31px;font-weight:650;letter-spacing:-.02em;margin:3px 0 1px}
.tile .sub{display:block;font-size:12px;color:var(--muted);margin:0}
.tile.crit .fig{color:var(--critical)}
.tile.cat  .fig{color:var(--catalogue)}
.tile.mod  .fig{color:var(--model)}

/* Antibiogram — status colour AND the letter, never colour alone */
.antibiogram{border-collapse:separate;border-spacing:3px;width:100%}
.antibiogram th{font-size:11px;text-transform:uppercase;letter-spacing:.05em;
                color:var(--muted);font-weight:600;padding:2px 4px}
.antibiogram th.sid{text-align:right;font-size:12.5px;color:var(--ink2);
                    text-transform:none;letter-spacing:0;white-space:nowrap}
/* The letter is the primary channel and the fill reinforces it, so the ink is
   chosen per fill rather than fixed: white on the critical red measures 4.76:1
   and near-black on the good green 6.27:1, where white on that green would be
   3.35:1 and fail at this size. */
td.ab{text-align:center;border-radius:6px;font-weight:700;font-size:14px;
      padding:9px 0;border-bottom:3px solid transparent}
td.ab.r{background:var(--critical);color:#fff}
td.ab.s{background:var(--good);color:#08240b}
td.ab.none{background:var(--grid);color:var(--muted)}
td.ab.l1{border-bottom-color:var(--catalogue)}
td.ab.l2{border-bottom-color:var(--model)}
td.ab span{cursor:help}
.legend{display:flex;gap:16px;flex-wrap:wrap;margin-top:12px;font-size:12.5px;color:var(--ink2)}
.legend i{display:inline-block;width:11px;height:11px;border-radius:3px;
          margin-right:6px;vertical-align:-1px;font-style:normal}
.legend .u{display:inline-block;width:16px;height:3px;border-radius:2px;
           margin-right:6px;vertical-align:3px}

/* per-drug cards */
.card{background:var(--surface);border:1px solid var(--border);border-radius:11px;
      padding:13px 16px;margin-bottom:9px}
.card>summary{cursor:pointer;list-style:none;display:flex;gap:12px;align-items:center;flex-wrap:wrap}
.card>summary::-webkit-details-marker{display:none}
.dr{font-weight:700;min-width:46px}
.call{font-weight:700;padding:1px 8px;border-radius:5px;color:#fff;font-size:12.5px}
.call.r{background:var(--critical)}.call.s{background:var(--good)}
.lyr{font-size:11.5px;padding:1px 8px;border-radius:99px;border:1px solid var(--border)}
.lyr.l1{color:var(--catalogue);border-color:var(--catalogue)}
.lyr.l2{color:var(--model);border-color:var(--model)}
.pr{font-size:12.5px;color:var(--ink2)}
.meter{display:inline-flex;align-items:center;gap:3px;margin-left:auto}
.meter .seg{width:11px;height:7px;border-radius:2px;background:var(--grid)}
.meter .seg.on{background:var(--ink2)}
.meter .mv{font-size:11.5px;color:var(--muted);margin-left:5px}
.meter .cap{font-size:11px;color:var(--warning)}

/* diverging contribution plot */
svg.wf{width:100%;height:auto;margin:12px 0 4px;overflow:visible}
.wfl{font-size:11.5px;fill:var(--ink2)}
.wfl.off{fill:var(--muted)}
.wfv{font-size:11px;fill:var(--muted)}
.wfa{font-size:10.5px;fill:var(--muted);text-transform:uppercase;letter-spacing:.05em}

table.graded,table.data{width:100%;border-collapse:collapse;font-size:13px;margin-top:9px}
table.graded th,table.data th{text-align:left;font-size:10.5px;text-transform:uppercase;
  letter-spacing:.06em;color:var(--muted);border-bottom:1px solid var(--grid);padding:5px 8px}
table.graded td,table.data td{padding:5px 8px;border-bottom:1px solid var(--grid)}
table.data td.r,table.graded td.r{color:var(--critical);font-weight:650}
table.data td.s,table.graded td.s{color:var(--good);font-weight:650}
.scroll{overflow-x:auto}

/* reasoning chain */
.chainwrap{margin-top:10px;border-top:1px solid var(--grid);padding-top:9px}
summary.mini{cursor:pointer;font-size:12.5px;color:var(--ink2);list-style:none}
summary.mini::-webkit-details-marker{display:none}
ol.chain{margin:10px 0 0;padding-left:20px;font-size:13px}
ol.chain li{margin-bottom:13px}
.stage{font-weight:650;text-transform:uppercase;font-size:11px;letter-spacing:.05em}
.kind{font-size:10.5px;color:var(--muted);margin-left:7px;border:1px solid var(--border);
      border-radius:99px;padding:0 6px}
ol.chain p{margin:4px 0}
.mech{color:var(--ink);border-left:3px solid var(--div_pos);padding-left:9px}
.assoc{color:var(--ink2);border-left:3px solid var(--muted);padding-left:9px}
.basis{font-size:11.5px;color:var(--muted)}
.note{border-left:3px solid var(--warning);padding:11px 14px;background:var(--surface);
      border-radius:8px;font-size:13px;color:var(--ink2);margin-bottom:14px}
.prov{font-size:12.5px;color:var(--ink2)}
.prov dt{font-weight:650;color:var(--ink);margin-top:7px;font-size:12px}
.prov dd{margin:1px 0 0 0}
footer{margin-top:34px;font-size:12px;color:var(--muted)}
@media print{
  body{background:#fff}
  .card,details{break-inside:avoid}
  details{open:true}
  .card>summary,summary.mini{list-style:none}
}
"""


def render_report(rows: list[dict], inputs: dict | None = None,
                  title: str = "mtb-resistotyper-ml report") -> str:
    """Assemble the whole result set into one offline HTML document."""
    inputs = inputs or {}
    by_sample: dict[str, list[dict]] = {}
    for r in rows:
        by_sample.setdefault(r["sample_id"], []).append(r)
    drugs = sorted({r["drug"] for r in rows})
    cat = next((r["provenance"].get("catalogue") for r in rows
                if r.get("provenance", {}).get("catalogue")), {}) or {}
    runner = next((r["provenance"].get("runner") for r in rows
                   if r.get("provenance", {}).get("runner")), "")
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    parts = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>{_e(title)}</title>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'>",
        "<style>",
        f":root{{{_vars(LIGHT)}}}",
        f"@media (prefers-color-scheme:dark){{:root:not([data-theme='light']){{{_vars(DARK)}"
        "color-scheme:dark;}}",
        f":root[data-theme='dark']{{{_vars(DARK)}color-scheme:dark;}}",
        CSS, "</style></head><body><div class='wrap'>",
        f"<h1>{_e(title)}</h1>",
        f"<p class='sub'>{len(by_sample)} isolate(s) · {len(rows)} drug calls · "
        f"generated {now}</p>",

        "<div class='note'><b>Research use only. Not a diagnostic.</b> Layer 2 "
        "models have not been prospectively validated and have no regulatory "
        "clearance. Do not use any value in this report to guide the treatment "
        "of a patient.</div>",

        _tiles(rows),

        "<h2>Susceptibility summary</h2>",
        "<div class='panel'>", _antibiogram(by_sample, drugs),
        "<div class='legend'>"
        f"<span><i style='background:{LIGHT['critical']}'></i>R — resistant</span>"
        f"<span><i style='background:{LIGHT['good']}'></i>S — susceptible</span>"
        f"<span><span class='u' style='background:var(--catalogue)'></span>"
        "graded by the WHO catalogue (layer 1)</span>"
        f"<span><span class='u' style='background:var(--model)'></span>"
        "answered by a model (layer 2)</span></div></div>",
    ]

    for sid, rs in by_sample.items():
        n1 = sum(1 for r in rs if r.get("layer") == 1)
        parts.append(f"<h2>{_e(sid)} — {n1} catalogue-graded, {len(rs) - n1} modelled</h2>")
        parts += [_drug_card(r) for r in sorted(rs, key=lambda x: x["drug"])]

    # The table view: the relief for every colour-carried distinction above, and
    # the thing a reader pastes into a spreadsheet.
    parts.append("<h2>All calls</h2><div class='panel scroll'><table class='data'>"
                 "<thead><tr><th>Sample</th><th>Drug</th><th>Call</th><th>Source</th>"
                 "<th>Layer</th><th>P(R)</th><th>Reliability</th><th>Tier</th>"
                 "<th>AUC</th><th>Reasoning</th></tr></thead><tbody>")
    for r in rows:
        rel, conf = r.get("reliability", {}), r.get("confidence", {})
        if r.get("layer") == 1:
            why = ";".join(f'{x["mutation"]}={x.get("call")}' for x in r.get("reasons", [])) or "—"
        else:
            why = ";".join(f'{x["mutation"]}({x["logit_contribution"]:+.2f})'
                           for x in r.get("reasons", [])) or "intercept only"
        parts.append(
            f"<tr><td>{_e(r['sample_id'])}</td><td>{_e(r['drug'])}</td>"
            f"<td class='{'r' if r['prediction'] == 'R' else 's'}'>{r['prediction']}</td>"
            f"<td>{_e(r.get('source', ''))}</td><td>{r.get('layer', '')}</td>"
            f"<td>{'' if r.get('p_resistant') is None else r['p_resistant']}</td>"
            f"<td>{rel.get('band', '')}/{rel.get('of', '')}</td>"
            f"<td>{_e(conf.get('tier') or '')}</td>"
            f"<td>{'' if conf.get('auc') is None else conf['auc']}</td>"
            f"<td>{_e(why)}</td></tr>")
    parts.append("</tbody></table></div>")

    parts.append("<h2>Provenance</h2><div class='panel prov'><dl>")
    parts.append(f"<dt>Runner</dt><dd>{_e(runner)}</dd>")
    if cat:
        parts.append(f"<dt>Layer 1 catalogue</dt><dd>{_e(cat.get('catalogue'))} "
                     f"version {_e(cat.get('version'))} · {_e(cat.get('grammar'))} · "
                     f"reference {_e(cat.get('reference'))}<br>"
                     f"<span class='basis'>file: {_e(cat.get('file'))}</span></dd>")
    card = next((r["provenance"].get("model_card") for r in rows
                 if r.get("provenance", {}).get("model_card")), {}) or {}
    if card:
        parts.append("<dt>Layer 2 models</dt><dd>" + " · ".join(
            f"{_e(k)}: {_e(v)}" for k, v in card.items()) + "</dd>")
    if inputs:
        src = inputs.get("source", {})
        parts.append("<dt>Input</dt><dd>" + _e(json.dumps(src)) + "</dd>")
        if inputs.get("lineage"):
            parts.append(f"<dt>Lineage covariate</dt><dd>{_e(inputs['lineage'])}</dd>")
    parts.append("</dl></div>")

    parts.append(
        "<footer>Layer 1 is the curated WHO catalogue read by piezo; Layer 2 is a "
        "versioned model bundle. They do not vote: a drug the catalogue can grade "
        "is answered by the catalogue and no model is consulted for it. "
        "Reproduce any row with <code>pip install mtb-resistotyper-ml</code>."
        "</footer></div></body></html>")
    return "".join(parts)
