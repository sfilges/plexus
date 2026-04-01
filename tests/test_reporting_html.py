"""Tests for HTML QC report generation (REPT-02)."""

from __future__ import annotations

import json

import pytest

from plexus.reporting.html_report import (
    generate_html_report,
    generate_html_report_from_data,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

MINIMAL_QC = {
    "tm_distribution": {
        "mean": 60.0,
        "std": 1.5,
        "min": 58.0,
        "max": 62.0,
        "per_primer": [
            {
                "junction": "GENE_A",
                "direction": "forward",
                "name": "GENE_A_fwd",
                "tm": 59.5,
            },
            {
                "junction": "GENE_A",
                "direction": "reverse",
                "name": "GENE_A_rev",
                "tm": 60.5,
            },
            {
                "junction": "GENE_B",
                "direction": "forward",
                "name": "GENE_B_fwd",
                "tm": 58.0,
            },
            {
                "junction": "GENE_B",
                "direction": "reverse",
                "name": "GENE_B_rev",
                "tm": 62.0,
            },
        ],
    },
    "sequence_flags": {
        "gc_high_threshold": 70.0,
        "gc_low_threshold": 30.0,
        "homopolymer_min_run": 4,
        "high_gc_count": 0,
        "low_gc_count": 0,
        "homopolymer_count": 1,
        "flagged_primers": [
            {
                "junction": "GENE_B",
                "direction": "forward",
                "name": "GENE_B_fwd",
                "gc": 35.0,
                "sequence": "AAAATGCATGCATGCATGC",
                "flags": ["homopolymer"],
            }
        ],
    },
    "cross_reactivity_matrix": {
        "dimer_threshold": 0.0,
        "matrix": {
            "GENE_A": {
                "GENE_B": {"min_dimer_score": -5.2, "interaction_count": 4},
            },
            "GENE_B": {
                "GENE_A": {"min_dimer_score": -5.2, "interaction_count": 4},
            },
        },
    },
    "primer_dimer_matrix": {
        "primer_labels": [
            "GENE_A_forward",
            "GENE_A_reverse",
            "GENE_B_forward",
            "GENE_B_reverse",
        ],
        "dimer_threshold": 0.0,
        "matrix": {
            "GENE_A_forward": {
                "GENE_A_reverse": -3.1,
                "GENE_B_forward": -4.5,
                "GENE_B_reverse": -5.2,
            },
            "GENE_A_reverse": {
                "GENE_A_forward": -3.1,
                "GENE_B_forward": -2.8,
                "GENE_B_reverse": -4.0,
            },
            "GENE_B_forward": {
                "GENE_A_forward": -4.5,
                "GENE_A_reverse": -2.8,
                "GENE_B_reverse": -1.5,
            },
            "GENE_B_reverse": {
                "GENE_A_forward": -5.2,
                "GENE_A_reverse": -4.0,
                "GENE_B_forward": -1.5,
            },
        },
    },
}

MINIMAL_SUMMARY = {
    "panel_name": "TestPanel",
    "genome": "hg38",
    "num_junctions": 3,
    "num_selected_pairs": 2,
    "best_multiplex_cost": 1234.56,
}


@pytest.fixture
def qc_output_dir(tmp_path):
    """Create a minimal pipeline output directory."""
    (tmp_path / "panel_qc.json").write_text(json.dumps(MINIMAL_QC))
    (tmp_path / "panel_summary.json").write_text(json.dumps(MINIMAL_SUMMARY))
    return tmp_path


@pytest.fixture
def full_output_dir(qc_output_dir):
    """Output dir with all optional files."""
    csv = (
        "Junction,Chrom,Junction_Start,Junction_End,Pair_ID,Forward_Seq,Reverse_Seq,"
        "Forward_Full_Seq,Reverse_Full_Seq,Forward_Tm,Reverse_Tm,Tm_Diff,"
        "Forward_Bound,Reverse_Bound,Forward_GC,Reverse_GC,Forward_Length,Reverse_Length,"
        "Forward_Genomic_Start,Forward_Genomic_End,Reverse_Genomic_Start,Reverse_Genomic_End,"
        "Amplicon_Length,Insert_Size,Pair_Penalty,Dimer_Score,Off_Target_Count,"
        "Specificity_Checked,On_Target_Detected,SNP_Count,SNP_Penalty,"
        "Forward_SNP_Count,Reverse_SNP_Count,"
        "Forward_Self_Any_Th,Reverse_Self_Any_Th,"
        "Forward_Self_End_Th,Reverse_Self_End_Th,"
        "Forward_Hairpin_Th,Reverse_Hairpin_Th,"
        "Forward_End_Stability,Reverse_End_Stability,"
        "Forward_Penalty,Reverse_Penalty,"
        "Amplicon_GC,Num_Candidate_Pairs,Cross_Dimer_Contribution\n"
        "GENE_A,chr1,100,100,GENE_A_fwd_rev,ATCG,GCTA,ATCG,GCTA,"
        "59.5,60.5,1.0,80,70,50.0,55.0,20,20,80,99,101,120,80,40,"
        "100.0,-1.5,0,True,True,0,0.0,0,0,"
        "30.0,25.0,20.0,18.0,15.0,12.0,3.5,3.2,0.5,0.6,50.0,5,1.5\n"
        "GENE_B,chr2,200,200,GENE_B_fwd_rev,TTTT,CCCC,TTTT,CCCC,"
        "58.0,62.0,4.0,75,65,35.0,65.0,22,22,180,201,201,222,90,46,"
        "120.0,-2.0,1,True,True,0,0.0,0,0,"
        "32.0,28.0,22.0,20.0,16.0,14.0,3.8,3.4,0.7,0.8,55.0,3,2.0\n"
    )
    (qc_output_dir / "selected_multiplex.csv").write_text(csv)

    ot_csv = (
        "Pair_ID,Junction,OT_Chrom,OT_F_Primer,OT_R_Primer,OT_Product_Size,OT_F_Start,OT_R_Start\n"
        "GENE_B_fwd_rev,GENE_B,chr5,SEQ_1,SEQ_2,150,5000,5150\n"
    )
    (qc_output_dir / "off_targets.csv").write_text(ot_csv)

    fj_csv = (
        "Junction,Chrom,Start,End,Error\nGENE_C,chr3,300,300,no valid primer pairs\n"
    )
    (qc_output_dir / "failed_junctions.csv").write_text(fj_csv)

    return qc_output_dir


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestGenerateHtmlReport:
    def test_generates_html_file(self, qc_output_dir):
        path = generate_html_report(qc_output_dir)
        assert path.exists()
        assert path.name == "panel_report.html"
        assert path.stat().st_size > 1000

    def test_html_contains_panel_name(self, qc_output_dir):
        path = generate_html_report(qc_output_dir, panel_name="MyTestPanel")
        html = path.read_text()
        # panel_name from summary takes precedence
        assert "TestPanel" in html

    def test_html_contains_plotly(self, qc_output_dir):
        path = generate_html_report(qc_output_dir)
        html = path.read_text()
        assert "plotly.js" in html.lower() or "Plotly.newPlot" in html

    def test_html_contains_tm_data(self, qc_output_dir):
        path = generate_html_report(qc_output_dir)
        html = path.read_text()
        assert "59.5" in html
        assert "GENE_A" in html

    def test_html_contains_heatmap_data(self, qc_output_dir):
        path = generate_html_report(qc_output_dir)
        html = path.read_text()
        assert "-5.2" in html
        assert "cross_reactivity_matrix" in html or "heatmap" in html.lower()

    def test_html_contains_primer_heatmap(self, qc_output_dir):
        path = generate_html_report(qc_output_dir)
        html = path.read_text()
        assert "primer-heatmap-chart" in html
        assert "Primer-Level Dimer Heatmap" in html

    def test_html_handles_missing_optional_files(self, tmp_path):
        (tmp_path / "panel_qc.json").write_text(json.dumps(MINIMAL_QC))
        path = generate_html_report(tmp_path, panel_name="Bare")
        html = path.read_text()
        assert "Bare" in html
        assert path.exists()

    def test_html_handles_empty_panel(self, tmp_path):
        empty_qc = {
            "tm_distribution": {
                "mean": None,
                "std": None,
                "min": None,
                "max": None,
                "per_primer": [],
            },
            "sequence_flags": {
                "gc_high_threshold": 70.0,
                "gc_low_threshold": 30.0,
                "homopolymer_min_run": 4,
                "high_gc_count": 0,
                "low_gc_count": 0,
                "homopolymer_count": 0,
                "flagged_primers": [],
            },
            "cross_reactivity_matrix": {"dimer_threshold": 0.0, "matrix": {}},
        }
        (tmp_path / "panel_qc.json").write_text(json.dumps(empty_qc))
        path = generate_html_report(tmp_path, panel_name="Empty")
        assert path.exists()
        html = path.read_text()
        assert "Empty" in html

    def test_raises_without_qc_json(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="panel_qc.json"):
            generate_html_report(tmp_path)

    def test_full_report_with_all_files(self, full_output_dir):
        path = generate_html_report(full_output_dir)
        html = path.read_text()
        # Amplicon chart should be present
        assert "amplicon-chart" in html
        # Off-targets table
        assert "GENE_B" in html
        assert "chr5" in html
        # Failed junctions table
        assert "GENE_C" in html
        assert "no valid primer pairs" in html
        # GC chart
        assert "gc-chart" in html


class TestGenerateHtmlReportFromData:
    def test_from_data_creates_file(self, tmp_path):
        output = tmp_path / "report.html"
        path = generate_html_report_from_data(
            qc_data=MINIMAL_QC,
            summary_data=MINIMAL_SUMMARY,
            selected_pairs_csv=None,
            off_targets_csv=None,
            failed_junctions_csv=None,
            output_path=output,
            panel_name="DataTest",
        )
        assert path.exists()
        html = path.read_text()
        assert "TestPanel" in html  # from summary_data

    def test_from_data_with_csv_strings(self, tmp_path):
        csv_str = (
            "Junction,Forward_GC,Reverse_GC,Amplicon_Length\nGENE_A,50.0,55.0,80\n"
        )
        output = tmp_path / "report.html"
        path = generate_html_report_from_data(
            qc_data=MINIMAL_QC,
            summary_data=None,
            selected_pairs_csv=csv_str,
            off_targets_csv=None,
            failed_junctions_csv=None,
            output_path=output,
            panel_name="CSVTest",
        )
        assert path.exists()
        html = path.read_text()
        assert "CSVTest" in html
