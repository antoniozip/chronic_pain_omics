"""Which metabolite names get resolved to ChEBI, and why singletons must.

The map builder originally resolved only names measured in two or more
studies, reasoning that a name seen once cannot reach a pooled estimate. It
can: the join is between accessions, not spellings. Two studies naming one
compound differently meet only after both names resolve, and under the old
filter neither name qualified. These tests pin the corrected selection.
"""
from __future__ import annotations


def test_default_resolves_every_name(metabolite_map_builder):
    """A singleton name is selected, because it may share an accession."""
    studies = {
        "erucamide": {"MTBLS1", "MTBLS2"},
        "13zdocosenamide": {"MTBLS3"},
    }
    selected = metabolite_map_builder.select_names(
        studies, min_studies=2, shared_only=False)
    assert set(selected) == {"erucamide", "13zdocosenamide"}


def test_shared_only_drops_singletons(metabolite_map_builder):
    """The old behaviour stays reachable, for reproducing the 502-row map."""
    studies = {
        "erucamide": {"MTBLS1", "MTBLS2"},
        "13zdocosenamide": {"MTBLS3"},
    }
    selected = metabolite_map_builder.select_names(
        studies, min_studies=2, shared_only=True)
    assert set(selected) == {"erucamide"}


def test_shared_only_honours_min_studies(metabolite_map_builder):
    studies = {
        "taurine": {"MTBLS1", "MTBLS2", "MTBLS3"},
        "creatine": {"MTBLS1", "MTBLS2"},
    }
    selected = metabolite_map_builder.select_names(
        studies, min_studies=3, shared_only=True)
    assert set(selected) == {"taurine"}


def test_selection_does_not_alias_the_input(metabolite_map_builder):
    """Mutating the result must not corrupt the caller's study index."""
    studies = {"taurine": {"MTBLS1"}}
    selected = metabolite_map_builder.select_names(
        studies, min_studies=2, shared_only=False)
    selected.pop("taurine")
    assert "taurine" in studies


def test_study_sets_are_preserved(metabolite_map_builder):
    """n_studies is written from these sets, so they must survive intact."""
    studies = {"creatine": {"MTBLS2774", "MTBLS5667"}}
    for shared_only in (True, False):
        selected = metabolite_map_builder.select_names(
            studies, min_studies=2, shared_only=shared_only)
        assert selected["creatine"] == {"MTBLS2774", "MTBLS5667"}


def test_norm_keeps_stereo_prefixes(metabolite_map_builder):
    """L-malic and D-malic are different molecules; normalisation must not
    merge them. An earlier estimate that stripped these prefixes reported 170
    reachable metabolites where the truth was 162."""
    norm = metabolite_map_builder.norm
    assert norm("L-Malic acid") != norm("D-Malic acid")
    assert norm("L-Malic acid") == norm("l-malic  ACID")


def test_unnamed_pattern_rejects_feature_clusters(metabolite_map_builder):
    """Unannotated peak ids must never become NAME: keys that join studies."""
    unnamed = metabolite_map_builder.UNNAMED
    for junk in ("cluster_12", "Unknown 4", "peak-7", "m123t45",
                 "metabolite_9", "compound 3", "417.22"):
        assert unnamed.match(junk), junk
    for real in ("Taurine", "L-Malic acid", "13Z-Docosenamide"):
        assert not unnamed.match(real), real


class TestNameVariants:
    """MAF rows do not always hold a compound name.

    Some hold a pipe-separated synonym blob. Queried whole it matches nothing
    and cost a 120s timeout apiece, which is what made the first full-recall
    run unviable; its last field is usually the name that resolves.
    """

    def test_plain_name_is_its_own_only_variant(self, metabolite_map_builder):
        assert metabolite_map_builder.name_variants("Taurine") == ["Taurine"]

    def test_pipe_blob_is_split(self, metabolite_map_builder):
        blob = ("5-methoxy-2,2-dimethyl-10-propyl-2H,8H-benzodipyran-8-one|"
                "5-Methoxy-2,2-dimethyl-10-propyl-2H,8H-benzodipyran-8-one|GUT-70")
        variants = metabolite_map_builder.name_variants(blob)
        assert "GUT-70" in variants
        assert blob not in variants

    def test_shortest_variant_is_tried_first(self, metabolite_map_builder):
        """A trivial name resolves where a full IUPAC string does not."""
        variants = metabolite_map_builder.name_variants(
            "N-(carboxymethyl)-N-methylglycine|Sarcosine|MFCD00008131")
        assert variants[0] == "Sarcosine"

    def test_blank_and_duplicate_parts_are_dropped(self, metabolite_map_builder):
        assert metabolite_map_builder.name_variants("Taurine||Taurine| ") == ["Taurine"]

    def test_empty_name_yields_no_query(self, metabolite_map_builder):
        for empty in ("", "   ", "|", " | "):
            assert metabolite_map_builder.name_variants(empty) == [], repr(empty)

    def test_order_is_deterministic_among_equal_lengths(self, metabolite_map_builder):
        """The map must be reproducible, so ties cannot depend on set order."""
        blob = "bbb|aaa|ccc"
        first = metabolite_map_builder.name_variants(blob)
        assert first == metabolite_map_builder.name_variants(blob)
        assert first == ["bbb", "aaa", "ccc"]   # first-seen order preserved

    def test_overlong_variants_are_not_queried(self, metabolite_map_builder):
        """OLS4 times out on full IUPAC strings rather than returning no match,
        so each costs the whole timeout and cannot succeed. ~50 such names
        accounted for most of a 2h17m run. The cap is set from the data: the
        longest name ever resolved through OLS4 is 43 characters."""
        mod = metabolite_map_builder
        iupac = ("5-Methoxy-2,2-dimethyl-6-<(E)-2-methylbut-2-enoyl>-10-propyl-"
                 "2H,8H-benzo<1,2-b:3,4-b'>dipyran-8-one")
        assert len(iupac) > mod.MAX_QUERY_LEN
        assert mod.name_variants(iupac + "|GUT-70") == ["GUT-70"]

    def test_a_name_that_is_only_overlong_yields_no_query(self, metabolite_map_builder):
        mod = metabolite_map_builder
        assert mod.name_variants("x" * (mod.MAX_QUERY_LEN + 1)) == []

    def test_a_name_at_the_limit_is_still_queried(self, metabolite_map_builder):
        mod = metabolite_map_builder
        name = "y" * mod.MAX_QUERY_LEN
        assert mod.name_variants(name) == [name]


class TestCacheOnlyMode:
    """OLS4 rate-limits a sustained client to roughly one lookup per 45s after
    a few thousand requests. Because the map is written only after the whole
    loop, a stalled run strands every name it has already resolved. Cache-only
    mode finishes such a run from what is already known."""

    def test_cached_hit_is_returned_without_network(self, metabolite_map_builder, monkeypatch):
        mod = metabolite_map_builder
        monkeypatch.setattr(mod, "CACHE_ONLY", True)
        monkeypatch.setattr(mod, "get_json", _explode)
        cache = {"Taurine": ["CHEBI:15891", "taurine", "ols4_label"]}
        assert mod.ols4_lookup("Taurine", cache) == ("CHEBI:15891", "taurine", "ols4_label")

    def test_uncached_name_is_unresolved_not_an_error(self, metabolite_map_builder, monkeypatch):
        """An unattempted name is recorded as unresolved, exactly as a failed
        lookup is, so it can be picked up by a later run."""
        mod = metabolite_map_builder
        monkeypatch.setattr(mod, "CACHE_ONLY", True)
        monkeypatch.setattr(mod, "get_json", _explode)
        assert mod.ols4_lookup("Never-Looked-Up", {}) is None

    def test_label_backfill_makes_no_request(self, metabolite_map_builder, monkeypatch):
        mod = metabolite_map_builder
        monkeypatch.setattr(mod, "CACHE_ONLY", True)
        monkeypatch.setattr(mod, "get_json", _explode)
        assert mod.chebi_label("CHEBI:99999", {}) == ""

    def test_network_is_used_when_the_flag_is_off(self, metabolite_map_builder, monkeypatch):
        """The guard must gate on the flag, not disable lookups outright."""
        mod = metabolite_map_builder
        monkeypatch.setattr(mod, "CACHE_ONLY", False)
        calls = []

        def fake(url, **params):
            calls.append(params.get("q"))
            return {"response": {"docs": [{"obo_id": "CHEBI:1", "label": "x",
                                           "synonym": []}]}}
        monkeypatch.setattr(mod, "get_json", fake)
        assert mod.ols4_lookup("x", {}) == ("CHEBI:1", "x", "ols4_label")
        assert calls == ["x"]


def _explode(*a, **k):
    raise AssertionError("network request made in cache-only mode")
