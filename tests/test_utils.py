from cp_multiomics.utils import get_logger, set_seed


def test_set_seed_runs():
    set_seed(42)


def test_get_logger():
    logger = get_logger(__name__)
    assert logger is not None
