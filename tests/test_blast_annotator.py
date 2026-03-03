import pandas as pd
import pytest

from plexus.blast.annotator import BlastResultsAnnotator


@pytest.fixture
def sample_blast_df():
    """
    Create a sample BLAST dataframe matching -outfmt 6 with columns:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen
    """
    cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen".split()
    data = [
        # 1. Perfect 3' match (from_3prime=True, length=20, qend=20)
        ["P1_F", "chr7", 100.0, 20, 0, 0, 1, 20, 100, 119, 0.0, 40.0, "plus", 20],
        # 2. Internal match (from_3prime=False, qend=15, length=10)
        ["P1_F", "chr7", 100.0, 10, 0, 0, 6, 15, 200, 209, 0.0, 30.0, "plus", 20],
        # 3. Mismatching 3' match (pident=90)
        ["P1_R", "chr7", 90.0, 20, 2, 0, 1, 20, 300, 319, 0.0, 25.0, "minus", 20],
        # 4. 5'-anchored partial match: qend(15) < qlen(22), so from_3prime=False
        ["P2_F", "chr22", 100.0, 15, 0, 0, 1, 15, 500, 514, 0.01, 20.0, "plus", 22],
        # 5. 3'-end partial match: qend(22) == qlen(22), so from_3prime=True
        ["P2_R", "chr22", 100.0, 15, 0, 0, 8, 22, 600, 614, 0.01, 20.0, "minus", 22],
    ]
    return pd.DataFrame(data, columns=cols)


def test_build_annotation_dict(sample_blast_df):
    annotator = BlastResultsAnnotator(sample_blast_df)
    annotator.build_annotation_dict(length_threshold=15, evalue_threshold=1.0)

    assert "from_3prime" in annotator.annotations
    assert "predicted_bound" in annotator.annotations


def test_add_annotations(sample_blast_df):
    annotator = BlastResultsAnnotator(sample_blast_df)
    annotator.build_annotation_dict(length_threshold=15, evalue_threshold=1.0)
    annotator.add_annotations()

    df = annotator.blast_df

    # Hit 0: Perfect 3' match
    assert bool(df.iloc[0]["from_3prime"]) is True
    assert bool(df.iloc[0]["length_pass_3prime"]) is True
    assert bool(df.iloc[0]["predicted_bound"]) is True

    # Hit 1: Internal match (qend=15, qlen=20)
    assert bool(df.iloc[1]["from_3prime"]) is False
    assert bool(df.iloc[1]["predicted_bound"]) is False

    # Hit 2: 3' match with 2 mismatches (pident=90, mismatch=2) — within tolerance
    assert bool(df.iloc[2]["length_pass_3prime"]) is True
    assert bool(df.iloc[2]["predicted_bound"]) is True

    # Hit 3: 5'-anchored partial match (qend=15, qlen=22): from_3prime must be False
    assert bool(df.iloc[3]["from_3prime"]) is False
    assert bool(df.iloc[3]["predicted_bound"]) is False

    # Hit 4: 3'-end partial match (qend=22, qlen=22): from_3prime must be True
    assert bool(df.iloc[4]["from_3prime"]) is True
    assert (
        bool(df.iloc[4]["length_pass_3prime"]) is True
    )  # length=15 >= threshold=15, pident=100
    assert bool(df.iloc[4]["predicted_bound"]) is True


def test_get_predicted_bound(sample_blast_df):
    annotator = BlastResultsAnnotator(sample_blast_df)
    annotator.build_annotation_dict()
    annotator.add_annotations()

    bound_df = annotator.get_predicted_bound()
    assert len(bound_df) == 3  # Hits 0, 2, 4
    assert all(bound_df["predicted_bound"])


def test_summarise_by_primer(sample_blast_df):
    annotator = BlastResultsAnnotator(sample_blast_df)
    annotator.build_annotation_dict()
    annotator.add_annotations()
    annotator.summarise_by_primer()

    summary = annotator.blast_primer_df
    assert "primer_name" in summary.columns
    assert "total_alignments" in summary.columns

    # P1_F has 2 hits, both are predicted_bound with default thresholds
    p1f = summary[summary["primer_name"] == "P1_F"].iloc[0]
    assert p1f["total_alignments"] == 2
    assert p1f["predicted_bound"] == 1  # Hit 0 is bound, Hit 1 is not. Correct.


def test_from_3prime_annotation_semantics():
    """
    from_3prime must be True only when the alignment reaches the 3' end of the
    PRIMER (qend == qlen), not when qend == alignment length.

    Two scenarios the bug gets wrong:
      - 5'-anchored partial: qend == alignment_length but qend < qlen → False
      - 3'-end partial:      qend == qlen but qend != alignment_length → True
    """
    cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen".split()
    data = [
        # 5'-anchored partial: alignment spans only 5' portion of primer
        ["SEQ_0", "chr1", 100.0, 15, 0, 0, 1, 15, 100, 114, 0.01, 20.0, "plus", 22],
        # 3'-end partial: alignment spans only 3' portion of primer
        ["SEQ_1", "chr1", 100.0, 15, 0, 0, 8, 22, 200, 214, 0.01, 20.0, "plus", 22],
        # Full-length: qend == length == qlen (both implementations agree)
        ["SEQ_2", "chr1", 100.0, 22, 0, 0, 1, 22, 300, 321, 0.0, 44.0, "plus", 22],
    ]
    df = pd.DataFrame(data, columns=cols)
    annotator = BlastResultsAnnotator(df)
    annotator.build_annotation_dict(length_threshold=15, evalue_threshold=1.0)
    annotator.add_annotations()

    # 5'-anchored: primer cannot extend → must NOT be flagged
    assert bool(df.iloc[0]["from_3prime"]) is False
    assert bool(df.iloc[0]["predicted_bound"]) is False

    # 3'-end partial: primer CAN extend → MUST be flagged
    assert bool(df.iloc[1]["from_3prime"]) is True
    assert bool(df.iloc[1]["predicted_bound"]) is True

    # Full-length: from_3prime=True, predicted_bound=True
    assert bool(df.iloc[2]["from_3prime"]) is True
    assert bool(df.iloc[2]["predicted_bound"]) is True


def test_length_pass_rejects_high_mismatch():
    """
    A 3' hit with more mismatches than max_mismatches should fail
    length_pass_3prime even if the alignment is long enough.
    """
    cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen".split()
    data = [
        # 3' hit with 3 mismatches — should fail length_pass_3prime
        ["SEQ_0", "chr1", 85.0, 20, 3, 0, 1, 20, 100, 119, 50.0, 15.0, "plus", 20],
        # 3' hit with 2 mismatches — should pass length_pass_3prime
        ["SEQ_1", "chr1", 90.0, 20, 2, 0, 1, 20, 200, 219, 50.0, 15.0, "plus", 20],
    ]
    df = pd.DataFrame(data, columns=cols)
    annotator = BlastResultsAnnotator(df)
    annotator.build_annotation_dict(
        length_threshold=15, evalue_threshold=1.0, max_mismatches=2
    )
    annotator.add_annotations()

    # 3 mismatches: from_3prime=True but length_pass_3prime=False
    assert bool(df.iloc[0]["from_3prime"]) is True
    assert bool(df.iloc[0]["length_pass_3prime"]) is False

    # 2 mismatches: from_3prime=True and length_pass_3prime=True
    assert bool(df.iloc[1]["from_3prime"]) is True
    assert bool(df.iloc[1]["length_pass_3prime"]) is True
