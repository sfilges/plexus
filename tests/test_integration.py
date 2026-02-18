"""End-to-end integration tests for the plexus pipeline.

These tests run the full pipeline (load -> merge -> design -> SNP check ->
optimize -> output) against small fixture files committed to tests/data/.

Run only integration tests:  pytest -m integration -v
Skip integration tests:      pytest -m "not integration"
"""

from __future__ import annotations

import shutil

import pytest

from plexus.pipeline import run_pipeline

# Check if BLAST+ tools are available on PATH
_BLAST_AVAILABLE = all(
    shutil.which(tool) is not None
    for tool in ("blastn", "makeblastdb", "blast_formatter")
)

pytestmark = pytest.mark.integration


@pytest.fixture(scope="class")
def pipeline_result(fixture_csv, fixture_fasta, fixture_vcf, tmp_path_factory):
    """Run the full pipeline once and share the result across the test class."""
    tmp_dir = tmp_path_factory.mktemp("integration_output")
    result = run_pipeline(
        input_file=fixture_csv,
        fasta_file=fixture_fasta,
        output_dir=tmp_dir,
        panel_name="integration_test",
        design_method="simsen",
        run_blast=False,
        snp_vcf=fixture_vcf,
        skip_snpcheck=False,
        padding=200,
    )
    return result


class TestFullPipeline:
    """Tests that share a single pipeline execution."""

    def test_pipeline_succeeds(self, pipeline_result):
        assert pipeline_result.success is True
        assert pipeline_result.errors == []

    def test_junction_merging(self, pipeline_result):
        # 3 input junctions -> 2 after merge (two CLTCL1 on chr22 merge)
        assert pipeline_result.num_junctions == 2

        names = [j.name for j in pipeline_result.panel.junctions]
        merged = [n for n in names if n.count("CLTCL1") >= 2]
        assert len(merged) == 1, f"Expected merged CLTCL1 junction, got names: {names}"

    def test_design_regions_extracted(self, pipeline_result):
        for junction in pipeline_result.panel.junctions:
            assert junction.design_region is not None
            # Point mutation with padding=200: design_region = 401bp
            assert len(junction.design_region) == 401
            # Sequence should be valid DNA
            assert set(junction.design_region) <= {"A", "C", "G", "T", "N"}

    def test_primer_pairs_designed(self, pipeline_result):
        assert pipeline_result.num_primer_pairs > 0

        for junction in pipeline_result.panel.junctions:
            assert junction.primer_pairs is not None
            assert len(junction.primer_pairs) > 0

            for pair in junction.primer_pairs:
                assert set(pair.forward.seq) <= {"A", "C", "G", "T"}
                assert set(pair.reverse.seq) <= {"A", "C", "G", "T"}

    def test_snp_check_applied(self, pipeline_result):
        assert "snp_checked" in pipeline_result.steps_completed

        # Verify the SNP check populated snp_count on all primer pairs.
        # With small fixtures, high-AF variants may not overlap primer sites,
        # so we check the mechanism ran (field set to int) rather than
        # requiring specific overlap counts.
        all_pairs = [
            pair
            for j in pipeline_result.panel.junctions
            if j.primer_pairs
            for pair in j.primer_pairs
        ]
        assert len(all_pairs) > 0
        for pair in all_pairs:
            assert isinstance(pair.snp_count, int)
            assert isinstance(pair.snp_penalty, float)

    def test_multiplex_optimization(self, pipeline_result):
        # One selected pair per effective junction (2 junctions after merge)
        assert len(pipeline_result.selected_pairs) == 2

    def test_output_files_created(self, pipeline_result):
        out = pipeline_result.output_dir
        expected = [
            "candidate_pairs.csv",
            "selected_multiplex.csv",
            "top_panels.csv",
            "panel_summary.json",
        ]
        for filename in expected:
            path = out / filename
            assert path.exists(), f"Missing output file: {filename}"
            assert path.stat().st_size > 0, f"Empty output file: {filename}"


class TestPipelineSkipSnpcheck:
    """Separate pipeline run with SNP check disabled."""

    def test_pipeline_skip_snpcheck(self, fixture_csv, fixture_fasta, tmp_path):
        result = run_pipeline(
            input_file=fixture_csv,
            fasta_file=fixture_fasta,
            output_dir=tmp_path / "skip_snp_output",
            panel_name="skip_snp_test",
            design_method="simsen",
            run_blast=False,
            skip_snpcheck=True,
            padding=200,
        )
        assert result.success is True
        assert "snp_check_skipped" in result.steps_completed
        assert "snp_checked" not in result.steps_completed


@pytest.mark.skipif(not _BLAST_AVAILABLE, reason="BLAST+ tools not found on PATH")
class TestFullPipelineWithBlast:
    """Integration tests with BLAST specificity check enabled.

    Requires BLAST+ (blastn, makeblastdb, blast_formatter) on PATH.
    Skipped automatically when not available.
    """

    @pytest.fixture(scope="class")
    def blast_pipeline_result(self, fixture_csv, fixture_fasta, fixture_vcf, tmp_path_factory):
        tmp_dir = tmp_path_factory.mktemp("blast_integration")
        result = run_pipeline(
            input_file=fixture_csv,
            fasta_file=fixture_fasta,
            output_dir=tmp_dir,
            panel_name="blast_test",
            run_blast=True,
            snp_vcf=fixture_vcf,
            padding=200,
        )
        return result

    def test_blast_step_completed(self, blast_pipeline_result):
        assert "specificity_checked" in blast_pipeline_result.steps_completed

    def test_specificity_checked_flag(self, blast_pipeline_result):
        for junction in blast_pipeline_result.panel.junctions:
            for pair in junction.primer_pairs:
                assert pair.specificity_checked is True

    def test_blast_output_files(self, blast_pipeline_result):
        blast_dir = blast_pipeline_result.output_dir / "blast"
        assert blast_dir.exists()
        assert (blast_dir / "all_primers.fasta").exists()

    def test_selected_pairs_have_off_target_info(self, blast_pipeline_result):
        for pair in blast_pipeline_result.selected_pairs:
            assert isinstance(pair.off_target_products, list)
