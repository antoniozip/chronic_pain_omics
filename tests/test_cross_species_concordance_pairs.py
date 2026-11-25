"""ortholog.concordance.rodent_pairs: the pairing behind Fig. 2C."""

import pandas as pd

from cp_multiomics.ortholog.concordance import rodent_pairs


def _table(rows):
    return pd.DataFrame(rows, columns=["animal_feature_id", "human_feature_id",
                                       "animal_yi", "animal_padj"])


def test_pairs_join_through_the_human_ortholog_not_the_symbol():
    mouse = _table([("Abc1", "ABC1", 0.5, 0.01), ("Xyz", "XYZ", -0.2, 0.5)])
    rat = _table([("Abc1l", "ABC1", 0.3, 0.02), ("Xyz", "OTHER", 0.1, 0.5)])
    pairs = rodent_pairs(mouse, rat)
    assert list(pairs["human_feature_id"]) == ["ABC1"]
    assert pairs["concordant"].tolist() == [True]
    assert pairs["both_significant"].tolist() == [True]


def test_a_human_gene_reached_twice_is_dropped_not_multiplied():
    mouse = _table([("A1", "DUP", 0.5, 0.5), ("A2", "DUP", -0.5, 0.5),
                    ("B", "ONE", 0.1, 0.5)])
    rat = _table([("A", "DUP", 0.5, 0.5), ("B", "ONE", -0.1, 0.5)])
    pairs = rodent_pairs(mouse, rat)
    assert list(pairs["human_feature_id"]) == ["ONE"]
    assert pairs["concordant"].tolist() == [False]
