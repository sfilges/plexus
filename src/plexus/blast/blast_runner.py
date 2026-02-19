import os
import subprocess

import pandas as pd

from plexus.utils.env import check_executable


def _check_blast_tools() -> None:
    """Verify that the required BLAST+ executables are available on PATH."""
    missing = [
        tool
        for tool in ("blastn", "makeblastdb", "blast_formatter")
        if not check_executable(tool)
    ]
    if missing:
        raise RuntimeError(
            f"BLAST+ tool(s) not found on PATH: {', '.join(missing)}. "
            "Install NCBI BLAST+ with: conda install -c bioconda blast  "
            "or: sudo apt-get install ncbi-blast+"
        )


# https://github.com/JasonAHendry/multiply/blob/master/src/multiply/blast/runner.py
class BlastRunner:
    BLAST_COLS = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen"

    def __init__(self, input_fasta, reference_fasta):
        """
        Interface for running BLAST:
        - Creates BLAST database
        - Runs BLAST in archive mode
        - Reformats to a table
        - Converts table to a pandas data frame

        params
            input_fasta: str
                Path to a .fasta file containing query sequences
                for BLAST search. In the context of MULTIPLY,
                these are primer sequences.
            reference_fasta: str
                Reference sequence against which queries are BLASTed.
        """
        self.input_fasta = input_fasta
        self.reference_fasta = reference_fasta

    def create_database(self):
        """
        Check if a `blast` database has already been generated for
        the `reference_fasta`, if not, create one.

        """
        _check_blast_tools()

        # Define database name
        if not self.reference_fasta.endswith(
            ".fasta"
        ) and not self.reference_fasta.endswith(".fa"):
            # Warn but proceed? Or assume valid?
            pass

        # Simple heuristic for DB name
        self.db_path = self.reference_fasta.rsplit(".", 1)[0]

        # Check if database already exists
        db_suffixes = [".nhr", ".nin", ".nsq"]
        if all([os.path.isfile(f"{self.db_path}{suff}") for suff in db_suffixes]):
            print(f"BLAST database '{self.db_path}' already exists.")
            return self

        # Create database
        subprocess.run(
            [
                "makeblastdb",
                "-in",
                self.reference_fasta,
                "-dbtype",
                "nucl",
                "-parse_seqids",
                "-out",
                self.db_path,
            ],
            check=True,
        )

        return self

    def run(self, output_archive, word_size=7):
        """
        Run blast, writing a BLAST archive to `output_archive`

        Note that we run with output format 11 `-outfmt 11` to
        produce the archive; from this format you can convert to
        all other formats.

        """

        # Run
        subprocess.run(
            [
                "blastn",
                "-db",
                self.db_path,
                "-query",
                self.input_fasta,
                "-word_size",
                str(word_size),
                "-outfmt",
                "11",
                "-out",
                output_archive,
            ],
            check=True,
        )

        # Save as instance variable, for reformattings
        self.output_archive = output_archive

        return self

    def reformat_output_as_table(self, output_table):
        """
        Reformat the output from `self.run()` to a table form,
        e.g. `-outfmt 6`.

        """

        # Run
        subprocess.run(
            [
                "blast_formatter",
                "-archive",
                self.output_archive,
                "-outfmt",
                f"6 {self.BLAST_COLS}",
                "-out",
                output_table,
            ],
            check=True,
        )

        # Save
        self.output_table = output_table

        return self

    def _load_as_dataframe(self):
        """
        Load tabular BLAST output as a pandas dataframe

        """

        # Load as a dataframe
        self.blast_df = pd.read_csv(
            self.output_table, sep="\t", names=self.BLAST_COLS.split(" ")
        )

    def get_dataframe(self):
        """
        Return a dataframe of tabular BLAST results

        """

        self._load_as_dataframe()

        return self.blast_df
