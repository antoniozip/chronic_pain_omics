from cp_multiomics.genomics.sample_size import (
    detect_cohort,
    is_burden_study,
    parse_sample_size,
)


def test_parse_sums_cases_and_controls():
    text = "1,362 European ancestry cases, 44,047 European ancestry controls"
    assert parse_sample_size(text) == 45409


def test_parse_handles_single_group():
    assert parse_sample_size("2,149 European ancestry individuals") == 2149


def test_parse_multiple_comma_groups():
    # three cohorts pooled
    assert parse_sample_size("275 cases, 15,749 controls, 1,000 replication") == 17024


def test_parse_returns_minus_one_when_no_numbers():
    assert parse_sample_size("European ancestry cases and controls") == -1
    assert parse_sample_size("") == -1


def test_detect_cohort_ukb_variants():
    assert detect_cohort("Knee pain for three months (UKB data field 3773)") == "UKB"
    assert detect_cohort("23,456 UK Biobank participants") == "UKB"
    assert detect_cohort("cases from the uk biobank cohort") == "UKB"


def test_detect_cohort_finngen_and_23andme():
    assert detect_cohort("450,000 FinnGen participants") == "FinnGen"
    assert detect_cohort("cohort provided by 23andMe, Inc.") == "23andMe"


def test_detect_cohort_none():
    assert detect_cohort("300 European ancestry female cases, 203 controls") == ""


def test_is_burden_study():
    assert is_burden_study("ICD10 M25.56: Pain in knee (Gene-based burden)") is True
    assert is_burden_study("Knee pain (UKB data field 3773) (Gene-based)") is True
    assert is_burden_study("Chronic knee pain") is False
    assert is_burden_study("Fibromyalgia") is False


def test_parse_ignores_digits_embedded_in_alphanumeric_tokens():
    # A regression dropping the \b anchors from _INT_RE would pull 37 out of
    # GRCh37 and 19 out of COVID19, silently corrupting the sample-size sum.
    assert parse_sample_size("build GRCh37, COVID19 substudy") == -1


def test_parse_extracts_standalone_number_next_to_words():
    # A genuine standalone integer token must still be summed.
    assert parse_sample_size("45409 European ancestry individuals") == 45409


# ---------------------------------------------------------------------------
# N is the METAL SAMPLESIZE weight and the cohort tie-break, so a parse error
# is a 2x weight error on one study rather than a cosmetic one. The parser used
# to sum every integer in the string.
# ---------------------------------------------------------------------------

def test_parse_prefers_a_leading_total_over_a_parenthetical_breakdown():
    # 12,345 total, broken down as 6,789 + 5,556: the old sum was 24,690.
    text = "12,345 European ancestry individuals (6,789 cases and 5,556 controls)"
    assert parse_sample_size(text) == 12345


def test_parse_ignores_an_age_range_and_a_cohort_count():
    # "aged 40-69" and "from 3 cohorts" contributed 40, 69 and 3: 2,112.
    assert parse_sample_size("2,000 individuals aged 40-69 from 3 cohorts") == 2000


def test_parse_ignores_a_ukb_data_field_id():
    # A field id is not a sample size. Returning -1 makes harmonize_sumstats
    # refuse the study (recorded as harmonize_failed) instead of weighting it by
    # 3,773. The old test asserted 3773, i.e. blessed the defect.
    assert parse_sample_size("UKB data field 3773") == -1


def test_detect_cohort_precedence_is_deterministic():
    # A string naming two cohorts must return the documented first-match (UKB).
    assert detect_cohort("UK Biobank cases with FinnGen replication") == "UKB"
    assert detect_cohort("FinnGen discovery, no biobank overlap") == "FinnGen"


def test_detect_cohort_is_case_insensitive_for_all_keys():
    assert detect_cohort("recruited via FINNGEN") == "FinnGen"
    assert detect_cohort("cohort from 23ANDME inc") == "23andMe"
