"""mtb-resistotyper-ml: a runner for versioned M. tuberculosis resistance models.

The models this tool consumes are not shipped with it. They live in their own
repository with their own release cadence, and are joined to this runner only by
the artefact contract: ``feature_schema.json`` fixes the input, ``model_card.json``
carries the operating range. Either side can move without the other.
"""

from mtb_resistotyper_ml.explain import explain
from mtb_resistotyper_ml.score import ModelBundle, score
from mtb_resistotyper_ml.vectorize import fired, vectorize

__version__ = "0.1.0"
__all__ = ["ModelBundle", "score", "explain", "vectorize", "fired", "__version__"]
