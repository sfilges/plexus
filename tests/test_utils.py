# ================================================================================
# Tests for utility functions (GC clamp, check_kmer, etc.)
# ================================================================================

from plexus.utils.utils import check_gc_clamp, check_kmer, find_max_poly_gc


class TestCheckGcClamp:
    """Tests for the 3' GC clamp filter."""

    def test_zero_gc_at_3prime_fails(self):
        """last 5 = ATATA, GC=0, too weak."""
        assert check_gc_clamp("ATGCATATA") is False

    def test_one_gc_at_3prime_passes(self):
        """last 5 = ATATC, GC=1, pass."""
        assert check_gc_clamp("ATGCATATC") is True

    def test_two_gc_at_3prime_passes(self):
        """last 5 = AGATC, GC=2, pass."""
        assert check_gc_clamp("ATGCAGATC") is True

    def test_three_gc_at_3prime_passes(self):
        """last 5 = GGATC, GC=3, pass (boundary)."""
        assert check_gc_clamp("ATGCGGATC") is True

    def test_four_gc_at_3prime_fails(self):
        """last 5 = GGCTC, GC=4, too stable."""
        assert check_gc_clamp("ATGCGGCTC") is False

    def test_five_gc_at_3prime_fails(self):
        """last 5 = GCGCG, GC=5, too stable."""
        assert check_gc_clamp("ATGCGCGCG") is False

    def test_short_sequence(self):
        """Sequence shorter than window handled gracefully."""
        # "ATGC" has 2 G/C in 4 bases (entire seq used as window)
        assert check_gc_clamp("ATGC") is True
        # "AAAA" has 0 G/C
        assert check_gc_clamp("AAAA") is False

    def test_lowercase_handled(self):
        """Lowercase bases should be counted correctly."""
        assert check_gc_clamp("atgcatatc") is True  # 1 G/C


class TestCheckKmerGcClamp:
    """Test gc_clamp parameter integration in check_kmer."""

    def test_gc_clamp_disabled_by_default(self):
        """gc_clamp=0 means GC clamp is not applied."""
        # Sequence with 0 G/C at 3' end — would fail clamp, but filter is off
        seq = "GCGCGATATAATATA"  # GC ~40%, last 5 = ATATA, GC=0
        assert check_kmer(seq, gc_clamp=0) is True

    def test_gc_clamp_enabled_rejects_weak_3prime(self):
        """gc_clamp=1 rejects primers with 0 G/C at 3' end."""
        seq = "GCGCGATATAATATA"  # last 5 = ATATA, GC=0
        assert check_kmer(seq, gc_clamp=1) is False

    def test_gc_clamp_enabled_accepts_good_3prime(self):
        """gc_clamp=1 accepts primers with 1-3 G/C at 3' end."""
        seq = "GCGCGATATAATATC"  # last 5 = ATATC, GC=1
        assert check_kmer(seq, gc_clamp=1) is True


class TestFindMaxPolyGc:
    """Tests for the G/C-specific poly-run filter."""

    def test_four_g_rejected(self):
        """GGGG (4 > 3) should be rejected."""
        assert find_max_poly_gc("ATGGGGAT", 3) is True

    def test_four_c_rejected(self):
        """CCCC (4 > 3) should be rejected."""
        assert find_max_poly_gc("ATCCCCAT", 3) is True

    def test_three_g_passes(self):
        """GGG (3 <= 3) should pass."""
        assert find_max_poly_gc("ATGGGATC", 3) is False

    def test_three_c_passes(self):
        """CCC (3 <= 3) should pass."""
        assert find_max_poly_gc("ATCCCATG", 3) is False

    def test_four_a_passes(self):
        """AAAA is not G/C, should pass."""
        assert find_max_poly_gc("GCAAAAGC", 3) is False

    def test_four_t_passes(self):
        """TTTT is not G/C, should pass."""
        assert find_max_poly_gc("GCTTTTGC", 3) is False

    def test_mixed_gc_passes(self):
        """GGCC has G-run=2 and C-run=2, both <= 3."""
        assert find_max_poly_gc("ATGGCCAT", 3) is False

    def test_short_sequence_passes(self):
        """Sequence shorter than threshold always passes."""
        assert find_max_poly_gc("GGG", 3) is False

    def test_lowercase_handled(self):
        """Lowercase G/C should be detected."""
        assert find_max_poly_gc("atggggatc", 3) is True


class TestCheckKmerPolyGc:
    """Test max_poly_gc integration in check_kmer."""

    def test_gggg_rejected_by_default(self):
        """A primer with GGGG is rejected by the default max_poly_gc=3."""
        # GGGGCCGGATGTTCTGCTG — the real DCAF12L2 primer
        seq = "GGGGCCGGATGTTCTGCTG"
        assert check_kmer(seq, max_poly_X=5, max_poly_gc=3) is False

    def test_ggg_passes(self):
        """A primer with GGG passes max_poly_gc=3."""
        seq = "GGGACCGGATGTTCTGCTG"
        assert check_kmer(seq, max_poly_X=5, max_poly_gc=3) is True

    def test_poly_gc_can_be_relaxed(self):
        """Setting max_poly_gc=4 allows GGGG."""
        seq = "GGGGCCGGATGTTCTGCTG"
        assert check_kmer(seq, max_poly_X=5, max_poly_gc=4) is True
