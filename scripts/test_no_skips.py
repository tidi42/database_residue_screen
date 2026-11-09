"""CI must execute all BLAST-gated paths."""
import sys
import unittest
import warnings


def main() -> int:
    warnings.simplefilter("error", ResourceWarning)
    suite = unittest.defaultTestLoader.discover("tests")
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    print(f"Executed={result.testsRun}; skipped={len(result.skipped)}")
    return 0 if result.wasSuccessful() and not result.skipped else 1


# Process-pool workers re-import this script by path under the forkserver/spawn start methods (the Linux default
# since Python 3.14); the guard keeps them from re-running the suite.
if __name__ == "__main__":
    sys.exit(main())
