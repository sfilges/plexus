# ================================================================================
# Tests for Panel QC Report (REPT-01)
# ================================================================================

from unittest.mock import MagicMock

import pytest

from plexus.reporting.qc import generate_panel_qc


def _make_primer(name, seq, tm, gc):
    p = MagicMock()
    p.name = name
    p.seq = seq
    p.tm = tm
    p.gc = gc
    return p


def _make_pair(
    fwd_name, fwd_seq, fwd_tm, fwd_gc, rev_name, rev_seq, rev_tm, rev_gc, selected=True
):
    pair = MagicMock()
    pair.selected = selected
    pair.forward = _make_primer(fwd_name, fwd_seq, fwd_tm, fwd_gc)
    pair.reverse = _make_primer(rev_name, rev_seq, rev_tm, rev_gc)
    return pair


def _make_junction(name, pairs):
    j = MagicMock()
    j.name = name
    j.primer_pairs = pairs
    return j


# ---------------------------------------------------------------------------
# Tm distribution
# ---------------------------------------------------------------------------


def test_tm_distribution_basic():
    pair = _make_pair("F1", "ACGTACGT", 60.0, 50.0, "R1", "TGCATGCA", 62.0, 50.0)
    junc = _make_junction("J1", [pair])
    result = generate_panel_qc([junc])
    td = result["tm_distribution"]
    assert td["mean"] == pytest.approx(61.0, abs=0.01)
    assert td["min"] == 60.0
    assert td["max"] == 62.0
    assert td["std"] is not None


def test_tm_distribution_single_primer_pair():
    pair = _make_pair("F1", "ACGTACGT", 65.0, 50.0, "R1", "TGCATGCA", 65.0, 50.0)
    junc = _make_junction("J1", [pair])
    result = generate_panel_qc([junc])
    td = result["tm_distribution"]
    assert td["mean"] == 65.0
    assert td["min"] == 65.0
    assert td["max"] == 65.0


def test_tm_distribution_empty():
    result = generate_panel_qc([])
    td = result["tm_distribution"]
    assert td["mean"] is None
    assert td["std"] is None
    assert td["min"] is None
    assert td["max"] is None
    assert td["per_primer"] == []


def test_tm_distribution_per_primer_entries():
    pair = _make_pair("FWD", "ACGTACGT", 60.0, 50.0, "REV", "TGCATGCA", 62.0, 50.0)
    junc = _make_junction("KRAS", [pair])
    result = generate_panel_qc([junc])
    per_primer = result["tm_distribution"]["per_primer"]
    assert len(per_primer) == 2
    names = {e["name"] for e in per_primer}
    assert names == {"FWD", "REV"}
    directions = {e["direction"] for e in per_primer}
    assert directions == {"forward", "reverse"}


def test_tm_distribution_std_single_value():
    """std should be None when only one primer exists (one junction, one pair but std needs >=2)."""
    pair = _make_pair("F1", "ACGT", 60.0, 50.0, "R1", "TGCA", 60.0, 50.0)
    junc = _make_junction("J1", [pair])
    result = generate_panel_qc([junc])
    # Two primers (forward + reverse) means stdev is computable
    assert result["tm_distribution"]["std"] is not None


# ---------------------------------------------------------------------------
# Sequence flags
# ---------------------------------------------------------------------------


def test_high_gc_flagged():
    pair = _make_pair("F1", "GCGCGCGC", 60.0, 75.0, "R1", "ACGTACGT", 60.0, 50.0)
    junc = _make_junction("J1", [pair])
    result = generate_panel_qc([junc])
    sf = result["sequence_flags"]
    assert sf["high_gc_count"] == 1
    assert sf["low_gc_count"] == 0
    flagged = sf["flagged_primers"]
    assert len(flagged) == 1
    assert "high_gc" in flagged[0]["flags"]
    assert flagged[0]["name"] == "F1"


def test_low_gc_flagged():
    pair = _make_pair("F1", "ATATATATAT", 60.0, 25.0, "R1", "ACGTACGT", 60.0, 50.0)
    junc = _make_junction("J1", [pair])
    result = generate_panel_qc([junc])
    sf = result["sequence_flags"]
    assert sf["low_gc_count"] == 1
    flagged = sf["flagged_primers"]
    assert any("low_gc" in f["flags"] for f in flagged)


def test_homopolymer_flagged():
    # sequence with a run of 4 identical bases
    pair = _make_pair("F1", "ACGTAAAACGT", 60.0, 45.0, "R1", "TGCATGCA", 60.0, 50.0)
    junc = _make_junction("J1", [pair])
    result = generate_panel_qc([junc])
    sf = result["sequence_flags"]
    assert sf["homopolymer_count"] == 1
    flagged = sf["flagged_primers"]
    assert any("homopolymer" in f["flags"] for f in flagged)


def test_no_flags_clean_primers():
    pair = _make_pair("F1", "ACGTACGT", 60.0, 50.0, "R1", "TGCATGCA", 60.0, 50.0)
    junc = _make_junction("J1", [pair])
    result = generate_panel_qc([junc])
    sf = result["sequence_flags"]
    assert sf["high_gc_count"] == 0
    assert sf["low_gc_count"] == 0
    assert sf["homopolymer_count"] == 0
    assert sf["flagged_primers"] == []


def test_multiple_flags_on_same_primer():
    # Low GC and homopolymer together
    pair = _make_pair("F1", "AAAATATA", 60.0, 12.5, "R1", "TGCATGCA", 60.0, 50.0)
    junc = _make_junction("J1", [pair])
    result = generate_panel_qc([junc])
    sf = result["sequence_flags"]
    flagged = sf["flagged_primers"]
    assert len(flagged) == 1
    assert "low_gc" in flagged[0]["flags"]
    assert "homopolymer" in flagged[0]["flags"]


def test_sequence_flags_thresholds_in_output():
    result = generate_panel_qc([], gc_high_threshold=65.0, gc_low_threshold=35.0)
    sf = result["sequence_flags"]
    assert sf["gc_high_threshold"] == 65.0
    assert sf["gc_low_threshold"] == 35.0


# ---------------------------------------------------------------------------
# Cross-reactivity matrix
# ---------------------------------------------------------------------------


def test_matrix_symmetric():
    pair_a = _make_pair("FA", "ACGTACGT", 60.0, 50.0, "RA", "TGCATGCA", 60.0, 50.0)
    pair_b = _make_pair("FB", "CCCCGGGG", 60.0, 75.0, "RB", "GGGGCCCC", 60.0, 75.0)
    junc_a = _make_junction("JA", [pair_a])
    junc_b = _make_junction("JB", [pair_b])
    result = generate_panel_qc([junc_a, junc_b])
    matrix = result["cross_reactivity_matrix"]["matrix"]
    assert "JA" in matrix
    assert "JB" in matrix["JA"]
    assert "JB" in matrix
    assert "JA" in matrix["JB"]
    # Symmetric: same cell object
    assert matrix["JA"]["JB"] == matrix["JB"]["JA"]


def test_matrix_cell_keys():
    pair_a = _make_pair("FA", "ACGTACGT", 60.0, 50.0, "RA", "TGCATGCA", 60.0, 50.0)
    pair_b = _make_pair("FB", "CCCCGGGG", 60.0, 75.0, "RB", "GGGGCCCC", 60.0, 75.0)
    junc_a = _make_junction("JA", [pair_a])
    junc_b = _make_junction("JB", [pair_b])
    result = generate_panel_qc([junc_a, junc_b])
    cell = result["cross_reactivity_matrix"]["matrix"]["JA"]["JB"]
    assert "min_dimer_score" in cell
    assert "interaction_count" in cell


def test_matrix_empty_single_junction():
    """With only one junction, no pairs are compared and matrix should be empty."""
    pair = _make_pair("F1", "ACGTACGT", 60.0, 50.0, "R1", "TGCATGCA", 60.0, 50.0)
    junc = _make_junction("J1", [pair])
    result = generate_panel_qc([junc])
    matrix = result["cross_reactivity_matrix"]["matrix"]
    assert matrix == {}


def test_matrix_dimer_threshold_in_output():
    result = generate_panel_qc([], dimer_threshold=-2.0)
    assert result["cross_reactivity_matrix"]["dimer_threshold"] == -2.0


def test_only_selected_pairs_used():
    """Unselected pairs should not appear in QC."""
    selected = _make_pair(
        "F_sel", "ACGTACGT", 62.0, 50.0, "R_sel", "TGCATGCA", 62.0, 50.0, selected=True
    )
    unselected = _make_pair(
        "F_no", "GGGGGGGG", 70.0, 100.0, "R_no", "CCCCCCCC", 70.0, 100.0, selected=False
    )
    junc = _make_junction("J1", [selected, unselected])
    result = generate_panel_qc([junc])
    names = {e["name"] for e in result["tm_distribution"]["per_primer"]}
    assert "F_sel" in names
    assert "R_sel" in names
    assert "F_no" not in names
    assert "R_no" not in names


# ---------------------------------------------------------------------------
# Primer-level dimer matrix
# ---------------------------------------------------------------------------


def test_primer_dimer_matrix_basic():
    pair_a = _make_pair("FA", "ACGTACGT", 60.0, 50.0, "RA", "TGCATGCA", 60.0, 50.0)
    pair_b = _make_pair("FB", "CCCCGGGG", 60.0, 75.0, "RB", "GGGGCCCC", 60.0, 75.0)
    junc_a = _make_junction("JA", [pair_a])
    junc_b = _make_junction("JB", [pair_b])
    result = generate_panel_qc([junc_a, junc_b])
    pdm = result["primer_dimer_matrix"]
    assert len(pdm["primer_labels"]) == 4
    matrix = pdm["matrix"]
    # All 6 unique pairs should have scores (C(4,2) = 6)
    all_scores = []
    for la in pdm["primer_labels"]:
        for lb in pdm["primer_labels"]:
            if la != lb:
                assert lb in matrix.get(la, {}), f"Missing {la} vs {lb}"
                all_scores.append(matrix[la][lb])
    assert len(all_scores) == 12  # 4*3 entries (symmetric)


def test_primer_dimer_matrix_symmetric():
    pair_a = _make_pair("FA", "ACGTACGT", 60.0, 50.0, "RA", "TGCATGCA", 60.0, 50.0)
    pair_b = _make_pair("FB", "CCCCGGGG", 60.0, 75.0, "RB", "GGGGCCCC", 60.0, 75.0)
    junc_a = _make_junction("JA", [pair_a])
    junc_b = _make_junction("JB", [pair_b])
    result = generate_panel_qc([junc_a, junc_b])
    matrix = result["primer_dimer_matrix"]["matrix"]
    for la in matrix:
        for lb in matrix[la]:
            assert matrix[la][lb] == matrix[lb][la], f"Asymmetric: {la} vs {lb}"


def test_primer_dimer_matrix_labels_format():
    pair = _make_pair("F1", "ACGTACGT", 60.0, 50.0, "R1", "TGCATGCA", 60.0, 50.0)
    junc = _make_junction("KRAS", [pair])
    result = generate_panel_qc([junc])
    labels = result["primer_dimer_matrix"]["primer_labels"]
    assert labels == ["KRAS_forward", "KRAS_reverse"]


def test_primer_dimer_matrix_single_junction():
    """Single junction should still have F vs R score."""
    pair = _make_pair("F1", "ACGTACGT", 60.0, 50.0, "R1", "TGCATGCA", 60.0, 50.0)
    junc = _make_junction("J1", [pair])
    result = generate_panel_qc([junc])
    pdm = result["primer_dimer_matrix"]
    assert len(pdm["primer_labels"]) == 2
    matrix = pdm["matrix"]
    assert "J1_reverse" in matrix.get("J1_forward", {})
    assert "J1_forward" in matrix.get("J1_reverse", {})


def test_primer_dimer_matrix_empty():
    result = generate_panel_qc([])
    pdm = result["primer_dimer_matrix"]
    assert pdm["primer_labels"] == []
    assert pdm["matrix"] == {}


def test_primer_dimer_matrix_with_tails():
    """Scores should differ when tail sequences are provided."""
    pair_a = _make_pair(
        "FA", "ACGTACGTACGT", 60.0, 50.0, "RA", "TGCATGCATGCA", 60.0, 50.0
    )
    pair_b = _make_pair(
        "FB", "CCCCGGGGCCCC", 60.0, 75.0, "RB", "GGGGCCCCGGGG", 60.0, 75.0
    )
    junc_a = _make_junction("JA", [pair_a])
    junc_b = _make_junction("JB", [pair_b])

    result_bare = generate_panel_qc([junc_a, junc_b])
    result_tailed = generate_panel_qc(
        [junc_a, junc_b],
        forward_tail="GGACACTCTTTCCCTACACGAC",
        reverse_tail="GTGACTGGAGTTCAGACGTGT",
    )

    bare_matrix = result_bare["primer_dimer_matrix"]["matrix"]
    tailed_matrix = result_tailed["primer_dimer_matrix"]["matrix"]
    # At least one score should differ
    diffs = 0
    for la in bare_matrix:
        for lb in bare_matrix[la]:
            if bare_matrix[la][lb] != tailed_matrix[la][lb]:
                diffs += 1
    assert diffs > 0


def test_save_primer_dimer_matrix_csv(tmp_path):
    from plexus.reporting.qc import save_primer_dimer_matrix_csv

    pair_a = _make_pair("FA", "ACGTACGT", 60.0, 50.0, "RA", "TGCATGCA", 60.0, 50.0)
    pair_b = _make_pair("FB", "CCCCGGGG", 60.0, 75.0, "RB", "GGGGCCCC", 60.0, 75.0)
    junc_a = _make_junction("JA", [pair_a])
    junc_b = _make_junction("JB", [pair_b])
    result = generate_panel_qc([junc_a, junc_b])

    csv_path = str(tmp_path / "dimer_matrix.csv")
    save_primer_dimer_matrix_csv(result["primer_dimer_matrix"], csv_path)

    import csv

    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        rows = list(reader)

    # Header row + 4 data rows
    assert len(rows) == 5
    # Header: empty + 4 labels
    assert len(rows[0]) == 5
    assert rows[0][0] == ""
    # Diagonal should be empty
    for i in range(1, 5):
        assert rows[i][i] == ""
    # Off-diagonal should be numeric
    assert float(rows[1][2]) != 0 or rows[1][2] != ""
