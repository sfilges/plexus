# ================================================================================
# Environment and dependency verification utilities
# ================================================================================

from __future__ import annotations

import importlib.metadata
import shutil
import subprocess
from pathlib import Path


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


# ---------------------------------------------------------------------------
# Tool version capture
# ---------------------------------------------------------------------------


def get_tool_version(name: str) -> str | None:
    """Return version string for a system tool, or ``None`` if absent.

    Tries ``<name> --version`` first, then ``<name> -version``.
    Returns the first non-empty line of stdout/stderr.
    """
    if not check_executable(name):
        return None
    for flag in ("--version", "-version"):
        try:
            result = subprocess.run(
                [name, flag],
                capture_output=True,
                text=True,
                timeout=5,
            )
            output = (result.stdout or result.stderr or "").strip()
            if output:
                return output.splitlines()[0]
        except (subprocess.TimeoutExpired, OSError):
            continue
    return "unknown"


def get_tool_versions(tools: list[str]) -> dict[str, str | None]:
    """Return ``{tool_name: version_string}`` for each tool (``None`` if missing)."""
    return {t: get_tool_version(t) for t in tools}


def get_plexus_version() -> str:
    """Return the installed plexus version."""
    try:
        return importlib.metadata.version("plexus")
    except importlib.metadata.PackageNotFoundError:
        from plexus.version import __version__

        return __version__


def get_primer3_version() -> str | None:
    """Return the installed primer3-py version, or ``None``."""
    try:
        return importlib.metadata.version("primer3-py")
    except importlib.metadata.PackageNotFoundError:
        return None


class ComplianceError(RuntimeError):
    """Raised when the environment does not meet compliance manifest requirements."""


def load_compliance_manifest() -> dict:
    """Load the bundled compliance manifest JSON."""
    import json
    from importlib import resources

    text = (
        resources.files("plexus.data")
        .joinpath("compliance_manifest.json")
        .read_text(encoding="utf-8")
    )
    return json.loads(text)


def validate_environment(need_blast: bool = False, need_snp: bool = False) -> dict:
    """Check tool versions against the compliance manifest.

    Returns a structured verdict dict (for provenance). Raises
    ``ComplianceError`` if any required tool is missing or the wrong version.

    Verdict values per tool: ``"pass"`` | ``"fail"`` | ``"missing"`` | ``"unparseable"``
    """
    import re

    manifest = load_compliance_manifest()
    tools_spec = manifest["tools"]

    tools_to_check: list[str] = []
    if need_blast:
        tools_to_check.extend(["blastn", "makeblastdb", "blast_formatter"])
    if need_snp:
        tools_to_check.append("bcftools")

    verdicts: dict[str, dict] = {}
    failures: list[str] = []

    for tool_name in tools_to_check:
        if tool_name not in tools_spec:
            continue
        spec = tools_spec[tool_name]
        expected = spec["exact_version"]
        pattern = spec["version_regex"]

        version_str = get_tool_version(tool_name)

        if version_str is None:
            verdicts[tool_name] = {
                "expected": expected,
                "actual": None,
                "verdict": "missing",
            }
            failures.append(f"  - {tool_name}: not found on PATH")
        else:
            match = re.search(pattern, version_str)
            if not match:
                verdicts[tool_name] = {
                    "expected": expected,
                    "actual": version_str,
                    "verdict": "unparseable",
                }
                failures.append(
                    f"  - {tool_name}: could not parse version from {version_str!r}"
                )
            else:
                actual = match.group(1)
                if actual == expected:
                    verdicts[tool_name] = {
                        "expected": expected,
                        "actual": actual,
                        "verdict": "pass",
                    }
                else:
                    verdicts[tool_name] = {
                        "expected": expected,
                        "actual": actual,
                        "verdict": "fail",
                    }
                    failures.append(
                        f"  - {tool_name}: expected exactly {expected}, found {actual!r}"
                    )

    if failures:
        n = len(failures)
        msg = (
            f"Environment validation failed — {n} tool(s) do not match the compliance manifest:\n"
            + "\n".join(failures)
        )
        raise ComplianceError(msg)

    return {
        "manifest_version": manifest.get("version", "unknown"),
        **verdicts,
    }


def check_disk_space(path: str | Path, threshold_gb: float = 2.0) -> bool:
    """Verify that the target directory has sufficient free disk space.

    Args:
        path: Path to the directory to check.
        threshold_gb: Minimum free space required in GB (default: 2.0).

    Returns:
        True if space is above threshold, False otherwise.
    """
    from loguru import logger

    path = Path(path)
    # Use parent if the path doesn't exist yet
    check_path = path if path.exists() else path.parent

    # If parent also doesn't exist, we can't check yet (will be created)
    # but we can try the highest existing ancestor
    while not check_path.exists() and check_path.parent != check_path:
        check_path = check_path.parent

    try:
        usage = shutil.disk_usage(check_path)
        free_gb = usage.free / (1024**3)

        if free_gb < threshold_gb:
            logger.warning(
                f"Low disk space: {free_gb:.2f} GB free on {check_path}. "
                f"Threshold is {threshold_gb} GB. Pipeline may fail during output writing."
            )
            return False
        return True
    except (OSError, ValueError) as e:
        logger.debug(f"Could not check disk space on {check_path}: {e}")
        return True  # Proceed anyway if check fails
