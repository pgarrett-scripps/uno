import json
import os
import random
import struct
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import torch

from uno_rt.artifact import load_model, save_model, rebase_output
from uno_rt.dataset import Row, peptide_hash, read_data, split_data
from uno_rt.training import predict, train

ROOT = Path(__file__).resolve().parents[4]
MODELS = ROOT / "models"


def write_data(path, seed=1, count=30, units="minutes"):
    rng = random.Random(seed)
    seqs = ["".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=12)) for _ in range(count)]
    path.write_text("sequence\trt\n" + "".join(f"{s}\t{10 + i * .2}\n" for i, s in enumerate(seqs)))
    return seqs


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name)
        torch.set_num_threads(2)

    def tearDown(self):
        self.tmp.cleanup()

    def test_input_contract(self):
        path = self.path / "data.tsv"
        for text in ["PEPTIDEK\t1\n", "sequence\trt\nPEPTIDEK\tnan\n", "sequence\trt\nM[Ox]K\t1\n",
                     "sequence\trt\n" + "A" * 51 + "\t1\n", "sequence\trt\nPEPTIDEK\t-1\n",
                     "sequence\trt\trun_id\nPEPTIDEK\t1\ta\nPEPTIDER\t2\tb\n"]:
            path.write_text(text)
            with self.assertRaises(ValueError): read_data(path, "minutes")
        path.write_text("sequence\trt\nPEPTIDEK\t1\nPEPTIDEK\t3\n")
        self.assertEqual(read_data(path, "minutes"), [Row("PEPTIDEK", 2)])

    def test_peptide_split_and_run_holdout(self):
        rows = [Row(s, i, "a") for i, s in enumerate(["AAAA", "CCCC", "DDDD", "EEEE", "FFFF", "GGGG"])]
        rows += [Row("AAAA", 5, "b"), Row("CCCC", 8, "b")]
        tr, va = split_data(rows, holdout_run="b")
        self.assertEqual({r.sequence for r in va}, {"AAAA", "CCCC"})
        self.assertFalse({r.sequence for r in tr} & {r.sequence for r in va})
        a, b = split_data(rows, seed=9)
        self.assertEqual({r.sequence for r in b}, {r.sequence for r in split_data(list(reversed(rows)), seed=9)[1]})
        with self.assertRaises(ValueError): split_data(rows, validation=rows)

    def test_v1_v2_roundtrip_and_normalization(self):
        for name in ["rt_cnn_chrono_ch32.urtc", "rt_cnn_chrono_global_ch32.urtc", "rt_cnn_chrono_global_v3_ch32.urtc"]:
            model, info, _ = load_model(MODELS / name)
            seqs = ["PEPTIDEK", "ACDEFGHIK", "WWWWWWK"]
            old = predict(model, seqs, info)
            new = rebase_output(model, info, mean=-6., std=2.3)
            np.testing.assert_allclose(old, predict(model, seqs, new), atol=1e-4)
            out = self.path / name
            save_model(out, model, new, {"rt_units": "irt"})
            restored, ri, metadata = load_model(out)
            np.testing.assert_allclose(old, predict(restored, seqs, ri), atol=1e-4)
            self.assertEqual(metadata["rt_units"], "irt")
            with self.assertRaises(FileExistsError): save_model(out, model, new, {})

    def test_corrupt_model_rejected(self):
        original = (MODELS / "rt_cnn_chrono_global_v3_ch32.urtc").read_bytes()
        path = self.path / "bad.urtc"
        for data in [original[:80], original + b"extra", original[:4] + struct.pack("<I", 99) + original[8:]]:
            path.write_bytes(data)
            with self.assertRaises(ValueError): load_model(path)

    def test_train_export_and_lineage(self):
        data, out = self.path / "data.tsv", self.path / "custom.urtc"
        write_data(data)
        report = train(data, out, units="minutes", method="test C18", epochs=1, config={"device": "cpu", "threads": 2})
        model, info, meta = load_model(out)
        self.assertEqual(meta["rt_units"], "minutes")
        self.assertEqual(info["ch"], 32)
        self.assertTrue(np.isfinite(predict(model, ["PEPTIDEK"], info)).all())
        self.assertLessEqual(report["validation"]["mae"], report["baseline"]["mae"] + 1e-5)
        with self.assertRaisesRegex(ValueError, "starting model"):
            train(data, self.path / "leak.urtc", units="minutes", method="test", init=out, epochs=1)
        repeat = self.path / "repeat.urtc"
        again = train(data, repeat, units="minutes", method="test C18", epochs=1, config={"device": "cpu", "threads": 2})
        self.assertEqual(report, again)
        self.assertEqual(len(meta["seen_peptide_hashes"]), 30)

    @unittest.skipUnless(os.environ.get("UNO_TEST_BINARY"), "set UNO_TEST_BINARY to run CLI/Rust parity integration")
    def test_cli_training_evaluation_parity_and_independence(self):
        binary = str(Path(os.environ["UNO_TEST_BINARY"]).resolve())
        data, out = self.path / "data.tsv", self.path / "custom.urtc"
        write_data(data)
        advanced = self.path / "advanced.json"
        advanced.write_text(json.dumps({"device": "cpu", "threads": 2, "architecture": {"ch": 4, "emb": 3, "head": 4, "k": 3, "blocks": [1]}}))
        subprocess.run([binary, "rt", "train", "--data", str(data), "--out", str(out), "--units", "minutes", "--method", "test C18",
                        "--from-scratch", "--epochs", "1", "--config", str(advanced), "--python", os.sys.executable], check=True, capture_output=True, text=True)
        inspected = json.loads(subprocess.check_output([binary, "rt", "inspect", str(out)]))
        self.assertEqual(inspected["version"], 3)
        fine_config = self.path / "fine.json"
        fine_config.write_text(json.dumps({"device": "cpu", "threads": 2}))
        adapted = self.path / "adapted.urtc"
        subprocess.run([binary, "rt", "train", "--data", str(data), "--out", str(adapted), "--units", "minutes", "--method", "test C18",
                        "--epochs", "1", "--config", str(fine_config), "--python", os.sys.executable], check=True, capture_output=True, text=True)
        _, _, adapted_meta = load_model(adapted)
        from uno_rt.dataset import file_hash
        self.assertEqual(adapted_meta["parent_sha256"], file_hash(MODELS / "rt_cnn_chrono_global_v3_ch32.urtc"))
        held = self.path / "held.tsv"
        write_data(held, seed=999, count=10)
        predictions = self.path / "predictions.tsv"
        result = subprocess.run([binary, "rt", "evaluate", "--model", str(out), "--data", str(held), "--units", "minutes", "--out-predictions", str(predictions)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        model, info, _ = load_model(out)
        rows = read_data(held, "minutes")
        expected = predict(model, [r.sequence for r in rows], info)
        actual = np.loadtxt(predictions, delimiter="\t", skiprows=1, usecols=2)
        np.testing.assert_allclose(expected, actual, atol=1e-4, rtol=1e-5)
        leak = subprocess.run([binary, "rt", "evaluate", "--model", str(out), "--data", str(data), "--units", "minutes"], capture_output=True, text=True)
        self.assertNotEqual(leak.returncode, 0)
        self.assertIn("overlap", leak.stderr)
        # Old models with unknown output scales require separate calibration.
        for name in ["rt_cnn_chrono_ch32.urtc", "rt_cnn_chrono_global_ch32.urtc", "rt_cnn_chrono_global_v3_ch32.urtc"]:
            old_path = MODELS / name
            old, oi, _ = load_model(old_path)
            pred_file = self.path / (name + ".tsv")
            subprocess.run([binary, "rt", "evaluate", "--model", str(old_path), "--data", str(held), "--units", "minutes", "--calibration-data", str(data), "--out-predictions", str(pred_file)], check=True, capture_output=True)
            tr = read_data(data, "minutes")
            from uno_rt.training import align_on_training
            slope, intercept = align_on_training(np.array([r.rt for r in tr]), predict(old, [r.sequence for r in tr], oi))
            expected = predict(old, [r.sequence for r in rows], oi) * slope + intercept
            actual = np.loadtxt(pred_file, delimiter="\t", skiprows=1, usecols=2)
            np.testing.assert_allclose(expected, actual, atol=2e-4, rtol=1e-5)


if __name__ == "__main__":
    unittest.main()
