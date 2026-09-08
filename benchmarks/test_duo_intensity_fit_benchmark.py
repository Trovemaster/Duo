import tempfile
import unittest
from pathlib import Path

import duo_intensity_fit_benchmark as difb


REFERENCE = """\
<- 2.0 A1 0.0 A1 S 1558.4360 0.0000 1558.4360 3.52296604E-03 1.08811122E-09 9.14750691E-30 1.23300000E-29 1 1 0 0.0 0.0 1 0 0 0.0 0.0
<- 2.0 A1 0.0 A1 S 1572.4510 0.0000 1572.4510 3.01853492E-03 9.74993537E-10 8.10835680E-30 1.40600000E-29 1 1 0 1.0 1.0 1 0 0 0.0 0.0
<- 1.0 A1 1.0 A1 Q 3093.2029 1560.3054 1532.8975 3.92766539E-03 1.86152438E-09 6.11379549E-30 7.53200000E-33 1 2 0 1.0 1.0 1 1 0 1.0 1.0
"""


class IntensityFitBenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary_directory.name)
        self.reference = self.directory / "reference.int"
        self.actual = self.directory / "actual.int"
        self.reference.write_text(REFERENCE, encoding="utf-8")

    def tearDown(self):
        self.temporary_directory.cleanup()

    def compare(self, text: str, value_rtol: float = 1.0e-6):
        self.actual.write_text(text, encoding="utf-8")
        return difb.compare_files(
            self.reference, self.actual, value_rtol=value_rtol
        )

    def test_identical_file_passes(self):
        result, problems, comments = self.compare(REFERENCE)
        self.assertTrue(result.passed)
        self.assertEqual(result.matched, 3)
        self.assertFalse(problems)
        self.assertFalse(comments)

    def test_line_order_does_not_matter(self):
        reversed_text = "\n".join(reversed(REFERENCE.splitlines())) + "\n"
        result, _problems, _comments = self.compare(reversed_text)
        self.assertTrue(result.passed)
        self.assertEqual(result.matched, 3)

    def test_relative_intensity_tolerance(self):
        within = REFERENCE.replace("9.14750691E-30", "9.14751100E-30")
        result, _problems, _comments = self.compare(within)
        self.assertTrue(result.passed)

        outside = REFERENCE.replace("9.14750691E-30", "9.14760000E-30")
        result, _problems, _comments = self.compare(outside)
        self.assertFalse(result.passed)
        self.assertEqual(result.calculated_intensity.changed, 1)

    def test_assignment_change_is_only_a_comment(self):
        changed = REFERENCE.replace(
            "1 1 0 0.0 0.0 1 0 0 0.0 0.0",
            "9 8 7 0.0 0.0 6 5 4 0.0 0.0",
            1,
        )
        result, problems, comments = self.compare(changed)
        self.assertTrue(result.passed)
        self.assertFalse(problems)
        self.assertEqual(result.assignment_comments, 1)
        self.assertEqual(len(comments), 1)

    def test_missing_record_fails(self):
        result, _problems, _comments = self.compare(
            "\n".join(REFERENCE.splitlines()[:-1]) + "\n"
        )
        self.assertFalse(result.passed)
        self.assertEqual(result.missing, 1)

    def test_rigorous_label_change_fails(self):
        changed = REFERENCE.replace("<- 1.0 A1 1.0 A1 Q", "<- 1.0 A1 1.0 A1 R")
        result, _problems, _comments = self.compare(changed)
        self.assertFalse(result.passed)
        self.assertEqual(result.missing, 1)
        self.assertEqual(result.extra, 1)

    def test_reference_intensity_must_match_exactly(self):
        changed = REFERENCE.replace("1.23300000E-29", "1.23400000E-29")
        result, _problems, _comments = self.compare(changed)
        self.assertFalse(result.passed)
        self.assertEqual(result.reference_intensity.changed, 1)


if __name__ == "__main__":
    unittest.main()
