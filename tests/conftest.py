"""Shared pytest fixtures and module loader for numbered pipeline scripts."""
import importlib.util
import sys
from pathlib import Path

import pytest
import yaml

PIPELINE_DIR = Path(__file__).parent.parent / "pipeline"
SCRIPTS_DIR = Path(__file__).parent.parent / "scripts"
CONF_DIR = Path(__file__).parent.parent / "conf"


def _load_by_path(path: Path, module_name: str):
    if module_name in sys.modules:
        return sys.modules[module_name]
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = mod
    spec.loader.exec_module(mod)
    return mod


def load_pipeline_module(filename: str):
    """Load a pipeline script by filename (handles leading-digit names)."""
    return _load_by_path(PIPELINE_DIR / filename,
                         f"_pipeline_{filename.replace('.', '_').replace('-', '_')}")


def load_script_module(filename: str):
    """Load a helper script from scripts/ the same way."""
    return _load_by_path(SCRIPTS_DIR / filename,
                         f"_script_{filename.replace('.', '_').replace('-', '_')}")


@pytest.fixture(scope="session")
def metabolite_map_builder():
    return load_script_module("build_metabolite_id_map.py")


@pytest.fixture(scope="session")
def search_cfg() -> dict:
    with open(CONF_DIR / "search" / "chronic_pain.yaml") as f:
        return yaml.safe_load(f)


@pytest.fixture(scope="session")
def pipeline_01():
    return load_pipeline_module("01_search_literature.py")


@pytest.fixture(scope="session")
def pipeline_02():
    return load_pipeline_module("02_screen.py")


@pytest.fixture(scope="session")
def pipeline_03():
    return load_pipeline_module("03_ingest_omics.py")


@pytest.fixture(scope="session")
def pipeline_04():
    return load_pipeline_module("04_harmonize.py")


@pytest.fixture(scope="session")
def pipeline_07():
    return load_pipeline_module("07_cross_species.py")


@pytest.fixture(scope="session")
def pipeline_08():
    return load_pipeline_module("08_figures.py")
