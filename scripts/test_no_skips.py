"""CI must execute all BLAST-gated paths."""
import sys
import unittest
import warnings
warnings.simplefilter("error", ResourceWarning)
suite = unittest.defaultTestLoader.discover("tests")
result = unittest.TextTestRunner(verbosity=2).run(suite)
print(f"Executed={result.testsRun}; skipped={len(result.skipped)}")
sys.exit(0 if result.wasSuccessful() and not result.skipped else 1)
