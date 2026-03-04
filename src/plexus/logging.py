# ================================================================================
# Logging configuration using loguru
# ================================================================================

import sys
from datetime import datetime

from loguru import logger


class _ConsoleFilter:
    """Togglable filter for the stderr log handler.

    When disabled, INFO messages are suppressed on the console while
    Rich progress bars are active.  File logging is unaffected.
    """

    def __init__(self):
        self.enabled = True

    def __call__(self, record):
        return self.enabled


_console_filter = _ConsoleFilter()

# Remove default handler
logger.remove()

# Add console handler with colored output (filtered via _console_filter)
logger.add(
    sys.stderr,
    format=(
        "<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | "
        "<cyan>{name}</cyan>:<cyan>{function}</cyan> - <level>{message}</level>"
    ),
    level="INFO",
    colorize=True,
    filter=_console_filter,
)


def suppress_console_logging():
    """Disable the stderr log handler (used while progress bars are active)."""
    _console_filter.enabled = False


def restore_console_logging():
    """Re-enable the stderr log handler."""
    _console_filter.enabled = True


def configure_file_logging(log_dir: str = ".", debug: bool = False) -> str:
    """
    Configure file logging with a timestamped log file.

    Args:
        log_dir: Directory to write log files to. Defaults to current directory.
        debug: If True, write DEBUG+ messages. Otherwise write INFO+ only.

    Returns:
        Path to the log file.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_filename = f"{log_dir}/plexus_{timestamp}.log"

    logger.add(
        log_filename,
        format=(
            "{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | "
            "{name}:{function}:{line} - {message}"
        ),
        level="DEBUG" if debug else "INFO",
        rotation="10 MB",
    )

    return log_filename


# Re-export logger for convenient imports
__all__ = [
    "logger",
    "configure_file_logging",
    "suppress_console_logging",
    "restore_console_logging",
]
