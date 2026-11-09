"""UTC logging for runs and attempts, and the advisory run lock."""
import logging
import os
import time
from contextlib import contextmanager
from pathlib import Path

from .common import StudyError
from .config import redactor

try:
    import fcntl
except ImportError:  # non-POSIX platforms run without the advisory run lock
    fcntl = None

PACKAGE_LOGGER = "database_residue_screen"
ENGINE_LOGGER = "database_residue_screen.study.engine"


class UtcFormatter(logging.Formatter):
    """ISO 8601 UTC timestamps; secret values are replaced by <redacted> in every formatted line."""
    converter = time.gmtime

    def __init__(self, secret_values=()):
        super().__init__("%(asctime)s.%(msecs)03dZ %(levelname)s %(short_name)s: %(message)s", "%Y-%m-%dT%H:%M:%S")
        self.redact = redactor(secret_values)

    def format(self, record: logging.LogRecord) -> str:
        record.short_name = record.name.removeprefix(PACKAGE_LOGGER + ".")
        return self.redact(super().format(record))


@contextmanager
def run_logging(path: Path, secret_values=(), console=None):
    """Engine messages go to the run log (and to ``console`` if given)."""
    logger = logging.getLogger(ENGINE_LOGGER)
    handlers = [logging.FileHandler(path, encoding="utf-8")] + ([logging.StreamHandler(console)] if console is not None else [])
    previous = logger.level
    logger.setLevel(logging.INFO)
    for handler in handlers:
        handler.setFormatter(UtcFormatter(secret_values))
        logger.addHandler(handler)
    try:
        yield logger
    finally:
        for handler in handlers:
            logger.removeHandler(handler)
            handler.close()
        logger.setLevel(previous)


@contextmanager
def attempt_logging(path: Path, secret_values=()):
    """Every package log record (engine and stage loggers) also goes to the attempt log while the attempt runs."""
    logger = logging.getLogger(PACKAGE_LOGGER)
    handler = logging.FileHandler(path, encoding="utf-8")
    handler.setFormatter(UtcFormatter(secret_values))
    previous = logger.level
    if previous == logging.NOTSET or previous > logging.INFO:
        logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    try:
        yield
    finally:
        logger.removeHandler(handler)
        handler.close()
        logger.setLevel(previous)


class RunLock:
    """Exclusive, non-blocking advisory lock on RUN_DIR/run.lock for the duration of one invocation."""

    def __init__(self, path: Path):
        self.path = path
        self.handle = None

    def __enter__(self):
        self.handle = open(self.path, "a+", encoding="utf-8")
        if fcntl is not None:
            try:
                fcntl.flock(self.handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError:
                self.handle.close()
                raise StudyError(f"Run directory {self.path.parent} is in use by another run-study process") from None
        self.handle.seek(0)
        self.handle.truncate()
        self.handle.write(f"{os.getpid()}\n")
        self.handle.flush()
        return self

    def __exit__(self, *details):
        if fcntl is not None:
            fcntl.flock(self.handle.fileno(), fcntl.LOCK_UN)
        self.handle.close()


def lock_held(path: Path) -> bool:
    if fcntl is None or not path.is_file():
        return False
    with open(path, encoding="utf-8") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_SH | fcntl.LOCK_NB)
        except OSError:
            return True
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    return False
