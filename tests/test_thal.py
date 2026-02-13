# ================================================================================
# Tests for Thermodynamic alignment model
# ================================================================================

import pytest

from plexus.designer.thal import oligotm, seqtm, symmetry


class TestThal:
    def test_seqtm_short_sequence(self):
        """Test seqtm function with a very short sequence."""
        seq = "ATGC"
        result = seqtm(seq)
        # A 4bp sequence has a very low Tm (primer3-py also gives -51.8)
        assert result.Tm == pytest.approx(-51.8, abs=1.0)

    def test_seqtm_primer_sequence_at_53C(self):
        """Test seqtm function with a realistic primer sequence at 53°C annealing.

        At 53°C (below Tm of 58.1°C), expect high binding (~90%).
        """
        seq = "CAGTGGCTCTATTGAATTTCTGTG"
        result = seqtm(
            seq,
            dna_conc=50.0,
            salt_conc=50.0,
            divalent_conc=1.5,
            dntp_conc=0.6,
            annealing_temp=53.0,
        )

        # Tm should match primer3-py calculation
        assert result.Tm == pytest.approx(58.1, abs=0.5)
        # Amount bound should be high (~90%) since we're below Tm
        assert result.bound == pytest.approx(90.2, abs=2.0)

    def test_seqtm_primer_sequence_at_60C(self):
        """Test seqtm function with a realistic primer sequence at 60°C annealing.

        At 60°C (above Tm of 58.1°C), expect lower binding (~31%).
        """
        seq = "CAGTGGCTCTATTGAATTTCTGTG"
        result = seqtm(
            seq,
            dna_conc=50.0,
            salt_conc=50.0,
            divalent_conc=1.5,
            dntp_conc=0.6,
            annealing_temp=60.0,
        )

        # Tm should match primer3-py calculation
        assert result.Tm == pytest.approx(58.1, abs=0.5)
        # Amount bound should be lower (~31%) since we're above Tm
        assert result.bound == pytest.approx(31.1, abs=2.0)

    def test_seqtm_high_gc_primer(self):
        """Test seqtm with a high GC content primer."""
        seq = "GCGCGCGCGCGCGCGCGCGC"  # 100% GC, 20bp
        result = seqtm(seq)
        # High GC primers have higher Tm
        assert result.Tm > 70.0

    def test_oligotm_basic(self):
        """Test oligotm function directly."""
        seq = "ACGTACGTACGTACGTACGT"  # 20bp, 50% GC
        result = oligotm(seq)
        # Should return a reasonable Tm for a 20-mer
        assert 50.0 < result.Tm < 70.0
        assert 0.0 < result.bound < 100.0

    def test_symmetry(self):
        """Test symmetry detection for self-complementary sequences."""
        # Self-complementary sequences
        assert symmetry("ATAT") is True
        assert symmetry("GCGC") is True
        assert symmetry("ACGT") is True

        # Non-self-complementary sequences
        assert symmetry("AAAA") is False
        assert symmetry("ACGTA") is False  # Odd length
        assert symmetry("AATTGG") is False

    def test_annealing_temperature_effect(self):
        """Test that lower annealing temp gives higher bound fraction."""
        seq = "CAGTGGCTCTATTGAATTTCTGTG"

        result_low_temp = seqtm(seq, annealing_temp=50.0)
        result_high_temp = seqtm(seq, annealing_temp=70.0)

        # Lower annealing temp should give higher bound fraction
        assert result_low_temp.bound > result_high_temp.bound
        # Tm should be the same regardless of annealing temp
        assert result_low_temp.Tm == pytest.approx(result_high_temp.Tm, abs=0.1)
