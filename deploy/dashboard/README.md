# Layer-2 evaluation dashboard

Tests the claim the two-layer design rests on: that the models earn their place
by answering where the curated catalogue cannot.

## The cohort

CRyPTIC v3.4.0 isolates that WHO-UCN-GTB-PCI-2023.5 graded Unknown or Fail,
restricted to those carrying a measured MIC. 68,144 sample-drug pairs are
ungraded; 13,267 have a phenotype; 7,767 of those fall on the eight deployed
drugs, and 862 of them are phenotypically resistant.

This is the only cohort on which the design can be tested rather than described.
A model evaluated on isolates the catalogue already grades is being asked a
question that was never its job.

## The finding

At the threshold the models ship with, Layer 2 recovers **6 of 862**.

That reads as a failure and is not one. AUC on this cohort reaches 0.83 for RIF
and 0.79 for EMB — the ranking is informative; the cut is in the wrong place. The
models were fitted where 38.5% of isolates are resistant and deployed where 11%
are, so their base rate sits below 0.5 and almost nothing crosses it.
Recalibrating on this cohort, holding specificity at 0.90, recovers **137**.

AMI (AUC 0.55) and KAN (0.51) are the real negatives: no threshold rescues a
ranking at the no-information line, and for those two drugs the second layer
adds nothing here.

## What it says about the artefact contract

`model_card.json` declares an operating RANGE — an AUC and a tier — and no
operating POINT. A consumer therefore inherits 0.5, which is correct only where
prevalence matches the training cohort. The card should carry a threshold, the
prevalence it assumes, and the sensitivity and specificity that buys.

## Rebuilding

The evaluation is a batch job, baked into the image with the code that produced
it, so a screenshot traces to one cohort build and one scoring run. A dashboard
that recomputed per request would show a different number every time the models
moved underneath it.

```
scripts/build_cohort.py  <db> <models> <out>   # SQL -> ideal-input instances
scripts/score_cohort.py  <cohort> <models> <out>   # Layer 2, catalogue not re-run
scripts/analyse.py       <results>             # AUC + ROC sweep + recalibrated point
```
