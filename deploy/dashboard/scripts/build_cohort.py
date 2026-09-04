#!/usr/bin/env python3
"""Assemble the isolates the WHO catalogue could not grade but which have a phenotype.

This is the only cohort on which the two-layer design can actually be tested. A
model evaluated on isolates the catalogue already grades is being asked a
question that was never its job; the honest test is the set where Layer 1
returned U or F and a measured MIC exists to check the answer against.

Emits one JSON instance per (sample, drug) in the runner's ideal-input shape.
Mutations are restricted to each drug's feature schema at query time rather than
after: a sample carries hundreds of mutations, almost none of which the model has
a feature for, and pulling them all turns a minute into an hour.
"""
import json, sys
from pathlib import Path
import duckdb

DB = sys.argv[1] if len(sys.argv) > 1 else "/home/abhinavsharma/mtb-build/cryptic.duckdb"
MODELS = Path(sys.argv[2] if len(sys.argv) > 2 else "/opt/models")
OUT = Path(sys.argv[3] if len(sys.argv) > 3 else "/home/abhinavsharma/l2-eval")
OUT.mkdir(parents=True, exist_ok=True)

drugs = sorted(d.name for d in MODELS.iterdir() if (d / "feature_schema.json").exists())
c = duckdb.connect(DB, read_only=True)
print(f"drugs: {drugs}", flush=True)

manifest = {}
for drug in drugs:
    schema = json.loads((MODELS / drug / "feature_schema.json").read_text())
    wanted = {f"{f['gene']}@{f['mutation']}" for f in schema["features"]
              if f["kind"] == "mutation"}
    cohort = c.execute("""
        SELECT p.UNIQUEID AS sample_id, p.PREDICTION AS catalogue_call,
               ph.BINARY_PHENOTYPE AS phenotype, ph.MIC AS mic,
               any_value(g.LINEAGE) AS lineage
        FROM predictions p
        JOIN ukmyc_phenotypes ph
          ON ph.UNIQUEID = p.UNIQUEID AND ph.DRUG = p.DRUG
         AND ph.BINARY_PHENOTYPE IN ('R','S')
        LEFT JOIN genomes g ON g.UNIQUEID = p.UNIQUEID
        WHERE p.DRUG = ? AND p.PREDICTION IN ('U','F')
        GROUP BY 1,2,3,4
    """, [drug]).fetchdf()
    if cohort.empty:
        print(f"  {drug}: no unclassified isolates with a phenotype", flush=True)
        continue

    ids = cohort["sample_id"].tolist()
    muts = c.execute("""
        SELECT UNIQUEID, GENE, MUTATION
        FROM mutations
        WHERE UNIQUEID IN (SELECT unnest(?::VARCHAR[]))
          AND GENE || '@' || MUTATION IN (SELECT unnest(?::VARCHAR[]))
        GROUP BY 1,2,3
    """, [ids, sorted(wanted)]).fetchdf()
    by_sample = {}
    for uid, gene, mut in muts.itertuples(index=False):
        by_sample.setdefault(uid, []).append({"gene": gene, "mutation": mut})

    rows = []
    for r in cohort.itertuples(index=False):
        lin = (r.lineage or "unknown")
        lin = f"lineage{lin}" if lin and str(lin)[0].isdigit() else str(lin)
        rows.append({
            "sample_id": r.sample_id, "drug": drug,
            "catalogue_call": r.catalogue_call, "phenotype": r.phenotype,
            "mic": r.mic, "lineage": lin,
            "instance": {"sample_id": r.sample_id, "nomenclature": "GARC",
                         "variants": by_sample.get(r.sample_id, []),
                         "covariates": {"lineage": lin}},
        })
    (OUT / f"{drug}.json").write_text(json.dumps(rows))
    n_r = sum(1 for x in rows if x["phenotype"] == "R")
    carried = sum(1 for x in rows if x["instance"]["variants"])
    manifest[drug] = {"n": len(rows), "n_resistant": n_r,
                      "n_with_modelled_mutation": carried,
                      "n_features": len(wanted)}
    print(f"  {drug}: {len(rows):5,} isolates  {n_r:4,} R  "
          f"{carried:5,} carry a modelled mutation  ({len(wanted)} features)", flush=True)

(OUT / "cohort_manifest.json").write_text(json.dumps(manifest, indent=2))
print("DONE", flush=True)
