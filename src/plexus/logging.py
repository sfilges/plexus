# ================================================================================
# Logging configuration using loguru
# ================================================================================

import sys
from datetime import datetime

from loguru import logger

# Remove default handler
logger.remove()

# Add console handler with colored output
logger.add(
    sys.stderr,
    format=(
        "<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | "
        "<cyan>{name}</cyan>:<cyan>{function}</cyan> - <level>{message}</level>"
    ),
    level="INFO",
    colorize=True,
)


def configure_file_logging(log_dir: str = ".") -> str:
    """
    Configure file logging with a timestamped log file.

    Args:
        log_dir: Directory to write log files to. Defaults to current directory.

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
        level="DEBUG",
        rotation="10 MB",
    )

    return log_filename


# Re-export logger for convenient imports
__all__ = ["logger", "configure_file_logging"]
