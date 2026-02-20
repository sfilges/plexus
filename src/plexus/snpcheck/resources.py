# ================================================================================
# SNP resource constants and re-exports
#
# Provides shared constants and re-exports used by the SNP checking subsystem.
#
# Author: Stefan Filges (stefan.filges@pm.me)
# Copyright (c) 2025-2026 Stefan Filges
# ================================================================================

from __future__ import annotations

# Re-export get_cache_dir from the central resources module so that existing
# code and test mocks targeting ``plexus.snpcheck.resources.get_cache_dir``
# continue to work.
from plexus.resources import (
    DEFAULT_DATA_DIR as DEFAULT_CACHE_DIR,  # noqa: F401
)
from plexus.resources import (
    ENV_DATA_DIR,  # noqa: F401
    get_cache_dir,  # noqa: F401
)

ENV_SNP_VCF = "PLEXUS_SNP_VCF"
