from pathlib import Path

import pytest

from cp_multiomics.vocabulary import PainVocabulary, load_pain_vocabulary

REPO_ROOT = Path(__file__).parent.parent


def test_default_vocabulary_loads_from_conf():
    vocab = load_pain_vocabulary()
    assert isinstance(vocab, PainVocabulary)
    assert len(vocab.all_terms) == 119        # 4 core + 47 conditions + 68 mechanisms


def test_no_term_was_lost_in_extraction():
    v = load_pain_vocabulary()
    total = len(v.core) + len(v.conditions) + len(v.mechanisms) + len(v.context)
    assert total == 157                        # the original PAIN_KEYWORDS length


def test_core_pain_terms_present():
    terms = {t.lower() for t in load_pain_vocabulary().all_terms}
    for expected in ("chronic pain", "neuropathic pain", "inflammatory pain",
                     "nociplastic pain"):
        assert expected in terms


def test_conditions_absent_from_the_old_four_term_query_are_present():
    terms = {t.lower() for t in load_pain_vocabulary().all_terms}
    for expected in ("fibromyalgia", "migraine", "back pain",
                     "complex regional pain", "sciatica", "trigeminal neuralgia"):
        assert expected in terms


def test_generic_anatomy_is_context_and_never_pain_defining():
    """'brain' must not make a maize-leaf paper look like a pain study."""
    v = load_pain_vocabulary()
    all_lower = {t.lower() for t in v.all_terms}
    for generic in ("brain", "cortex", "hippocampus", "csf", "pons", "inflammation"):
        assert generic not in all_lower
    context_lower = {t.lower() for t in v.context}
    assert {"brain", "cortex", "csf", "pons", "inflammation"} <= context_lower


def test_matches_ignores_context_terms():
    v = load_pain_vocabulary()
    assert v.matches("Proteomic profiling of rat brain cortex") == ()


def test_all_terms_are_deduplicated_and_order_preserving():
    vocab = PainVocabulary(core=("a", "b"), conditions=("b", "c"), mechanisms=("a",))
    assert vocab.all_terms == ("a", "b", "c")


def test_all_terms_excludes_context():
    vocab = PainVocabulary(core=("chronic pain",), conditions=(), mechanisms=(),
                           context=("brain",))
    assert vocab.all_terms == ("chronic pain",)


def test_query_string_quotes_multiword_terms():
    vocab = PainVocabulary(core=("chronic pain",), conditions=("sciatica",), mechanisms=())
    assert vocab.query_string() == '"chronic pain" OR "sciatica"'
    assert vocab.query_string(quote=False) == "chronic pain OR sciatica"


def test_matches_is_case_insensitive():
    vocab = PainVocabulary(core=("chronic pain",), conditions=("Sciatica",), mechanisms=())
    assert vocab.matches("A study of CHRONIC PAIN in rats") == ("chronic pain",)
    assert vocab.matches("sciatica cohort") == ("Sciatica",)
    assert vocab.matches("maize leaf phosphoproteome") == ()


def test_vocabulary_is_immutable():
    vocab = load_pain_vocabulary()
    with pytest.raises(AttributeError):
        vocab.core = ()  # type: ignore[misc]


def test_vocabulary_defined_exactly_once_in_the_tree():
    """The 157-term list must not be re-typed in scripts/."""
    script = (REPO_ROOT / "scripts" / "analyze_proteomics.py").read_text()
    assert "PAIN_KEYWORDS = [" not in script
    assert "load_pain_vocabulary" in script


def test_short_acronyms_do_not_match_inside_longer_words():
    v = load_pain_vocabulary()
    assert v.matches("Proteomic analysis of occipital cortex synaptosomes") == ()
    assert v.matches("Venom gland proteome of Vipera ammodytes") == ()
    assert v.matches("CFAP53 knockout ciliary proteome") == ()
    assert v.matches("DRG2 GTPase complex in hepatocellular carcinoma") == ()


def test_short_acronyms_still_match_as_whole_words():
    v = load_pain_vocabulary()
    assert "DRG" in v.matches("DRG neurons after CCI")
    assert "CCI" in v.matches("DRG neurons after CCI")
    assert "DRG" in v.matches("ipsilateral DRGs were collected")   # plural preserved
    assert "SNI" in v.matches("SNI model of neuropathic pain")


def test_long_terms_still_match_as_substrings():
    v = load_pain_vocabulary()
    # 'nociceptive' (>4 chars) must still be found inside a larger token context
    assert "hyperalgesia" in v.matches("thermal hyperalgesia was measured")
    assert "chronic pain" in v.matches("RNA-seq of chronic pain")


def test_junk_records_still_match_nothing():
    v = load_pain_vocabulary()
    for junk in ("Direct evidence of milk consumption from ancient dental calculus",
                 "Surface Proteins of Listeria monocytogenes",
                 "maize leaf salt-responsive phosphoproteome"):
        assert v.matches(junk) == ()


def test_vocabulary_partition_invariants():
    v = load_pain_vocabulary()
    assert len(v.core) == 4
    assert len(v.conditions) == 47
    assert len(v.mechanisms) == 68
    assert len(v.context) == 38
    assert len(v.all_terms) == 119
    every = v.core + v.conditions + v.mechanisms + v.context
    assert len(every) == 157
    assert len(set(every)) == 157, "a term is duplicated across groups"
    assert all(t == t.strip() and t for t in every), "blank or padded term"
