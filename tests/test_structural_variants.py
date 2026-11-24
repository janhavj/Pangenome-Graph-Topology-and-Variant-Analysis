from pathlib import Path
import pandas as pd
import pytest

RESULTS = Path("results/tables/structural_variants.csv")


@pytest.fixture(scope="module")
def variants():
    assert RESULTS.exists(), (
        "Run python src/03_detect_structural_variants.py first"
    )
    return pd.read_csv(RESULTS).set_index("sample")


def test_three_variants_detected(variants):
    assert len(variants) == 3


def test_expected_samples(variants):
    assert set(variants.index) == {
        "insertion",
        "deletion",
        "inversion",
    }


@pytest.mark.parametrize(
    "sample,sv_type,start,end,length_change",
    [
        ("insertion", "INS", 30, 30, 10),
        ("deletion", "DEL", 20, 30, -10),
        ("inversion", "INV", 20, 30, 0),
    ],
)
def test_variant_properties(
    variants, sample, sv_type, start, end, length_change
):
    variant = variants.loc[sample]

    assert variant["variant_type"] == sv_type
    assert variant["ref_start"] == start
    assert variant["ref_end"] == end
    assert variant["length_change"] == length_change


def test_inversion_orientation(variants):
    variant = variants.loc["inversion"]

    assert variant["ref_nodes"] == "C+"
    assert variant["alt_nodes"] == "C-"


def test_insertion_segment(variants):
    assert variants.loc["insertion", "alt_nodes"] == "E+"


def test_deletion_segment(variants):
    assert variants.loc["deletion", "ref_nodes"] == "C+"
