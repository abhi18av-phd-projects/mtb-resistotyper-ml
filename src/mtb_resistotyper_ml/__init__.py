"""mtb-resistotyper-ml: a two-layer runner for M. tuberculosis resistance calls.

Layer 1 is the curated WHO catalogue, read by piezo. Layer 2 is a versioned model
bundle. They do not vote: a drug the catalogue can grade is answered by the
catalogue and no model is consulted for it.

The models are not shipped with this tool. They live in their own repository with
their own release cadence, joined to the runner only by the artefact contract:
``feature_schema.json`` fixes the input, ``model_card.json`` carries the operating
range. Either side can move without the other.
"""

from mtb_resistotyper_ml.version import __version__
from mtb_resistotyper_ml import catalogue
from mtb_resistotyper_ml.explain import explain
from mtb_resistotyper_ml.score import ModelBundle, reliability, score
from mtb_resistotyper_ml.vectorize import fired, vectorize
from mtb_resistotyper_ml.layers import catalogue_ready, resolve, resolve_all
from mtb_resistotyper_ml.ingest import (BcftoolsUnavailable, ReferenceUnavailable,
                                        instances_from_vcf, iter_instances_from_vcf)

__all__ = ["ModelBundle", "score", "reliability", "explain", "vectorize", "fired",
           "catalogue", "catalogue_ready", "resolve", "resolve_all",
           "instances_from_vcf", "iter_instances_from_vcf",
           "ReferenceUnavailable", "BcftoolsUnavailable", "__version__"]
