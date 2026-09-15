import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from matrix_koth070_study import Study


class DriverTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.study = Study(Path(self.tmp.name))

    def test_consensus_search_never_expands_members(self):
        matrix = self.study.out / 'aligned'
        matrix.mkdir()
        for name in ['matrix_manifest.json', 'feature_observations.tsv', 'lfq_matrix.long.tsv']:
            (matrix / name).write_text('fixture')
        with patch.object(self.study, 'job') as job:
            self.study.search(({'dataset': 'astral'}, 'consensus', matrix, 0, 'search'))
        _, argv, inputs, outputs, kind = job.call_args.args
        self.assertIn('--matrix-long', argv)
        self.assertNotIn('--matrix-joint', argv)
        self.assertNotIn('--features', argv)
        self.assertEqual(set(matrix.iterdir()), set(inputs) & set(matrix.iterdir()))
        self.assertEqual(len(outputs), 3)
        self.assertEqual(kind, 'search')

    def test_individual_search_uses_only_original_table(self):
        source = self.study.out / 'run/features.tsv'
        with patch.object(self.study, 'job') as job:
            self.study.search(({'dataset': 'astral'}, 'separate', source, 2, 'search'))
        _, argv, inputs, outputs, _ = job.call_args.args
        self.assertIn('--features', argv)
        self.assertNotIn('--matrix-long', argv)
        self.assertIn(source, inputs)
        self.assertEqual(len(outputs), 1)

    def test_resume_rejects_tampered_successful_output(self):
        lib = self.study.out / 'lib'
        lib.mkdir()
        (lib / 'lib_lightgbm.so').write_bytes(b'fingerprint fixture')
        output = self.study.out / 'result.txt'
        argv = [sys.executable, '-c', 'from pathlib import Path; import sys; Path(sys.argv[1]).write_text("ok")', output]
        self.study.job('tiny', argv, [], [output], 'test')
        self.study.job('tiny', argv, [], [output], 'test')
        self.assertEqual(len(list((self.study.out / 'jobs/tiny').glob('attempt-*'))), 1)
        output.write_text('changed')
        with self.assertRaises(AssertionError):
            self.study.job('tiny', argv, [], [output], 'test')
        record = json.loads((self.study.out / 'jobs/tiny/receipt.json').read_text())
        self.assertEqual(record['returncode'], 0)


if __name__ == '__main__':
    unittest.main()
