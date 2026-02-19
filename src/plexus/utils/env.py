# ================================================================================
# Environment and dependency verification utilities
# ================================================================================

import shutil


def check_executable(name: str) -> bool:
    """Check if an executable exists on the PATH."""
    return shutil.which(name) is not None


def get_missing_tools(need_blast: bool = False, need_snp: bool = False) -> list[str]:
    """
    Identify missing system dependencies based on run requirements.

    Args:
        need_blast: Whether BLAST+ tools are required.
        need_snp: Whether bcftools is required for SNP checking.

    Returns:
        List of missing executable names.
    """
    missing = []
    if need_blast:
        # We need both for a full run (db creation, searching, and formatting)
        for tool in ("blastn", "makeblastdb", "blast_formatter"):
            if not check_executable(tool):
                missing.append(tool)
    if need_snp:
        if not check_executable("bcftools"):
            missing.append("bcftools")
    return missing
