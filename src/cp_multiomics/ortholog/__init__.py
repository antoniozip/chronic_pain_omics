"""Ortholog mapping package.

Uses MyGene.info REST API to map gene symbols between species.
Results are cached locally to avoid repeat API calls.
"""

from .hcop import OrthologMapper

__all__ = ["OrthologMapper"]
