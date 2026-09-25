import unittest

from database_residue_screen.contracts import ContractError
from database_residue_screen.database_integrity import assert_single_build


class DatabaseVolumeTests(unittest.TestCase):
    def test_mixed_build_is_blocked(self):
        first = dict(metadata=dict(build_date="2025-06-04"))
        second = dict(metadata=dict(build_date="2026-05-29"))
        assert_single_build([first, first])
        for records in ([], [first, second]):
            with self.assertRaises(ContractError):
                assert_single_build(records)


if __name__ == "__main__":
    unittest.main()