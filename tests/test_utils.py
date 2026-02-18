# ================================================================================
# Tests for utility functions (GC clamp, check_kmer, etc.)
# ================================================================================

from plexus.utils.utils import check_gc_clamp, check_kmer


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
