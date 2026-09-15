import csv
import os
import random
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import torch
from uno_rt.artifact import load_model
from uno_rt.dataset import split_data
from uno_rt.im import read_data, train
from uno_rt.training import predict

class ImTests(unittest.TestCase):
    def test_charge_grouping_and_strict_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'data.tsv'
            p.write_text('sequence\tcharge\tim\nAAAA\t2\t1\nAAAA\t2\t3\nAAAA\t3\t4\nCCCC\t2\t1\nDDDD\t2\t2\nEEEE\t3\t2\n')
            rows=read_data(p,'inverse_mobility')
            self.assertEqual(len(rows),5)
            self.assertEqual(rows[0].rt,2)
            tr,va=split_data(rows)
            self.assertFalse({r.sequence for r in tr}&{r.sequence for r in va})
            for z,v in [('0','1'),('2','nan'),('2','-1'),('2.5','1')]:
                p.write_text(f'sequence\tcharge\tim\nAAAA\t{z}\t{v}\n')
                with self.assertRaises(ValueError):read_data(p,'inverse_mobility')

    def test_train_adapt_parity_and_lineage(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp); rng=random.Random(174)
            seqs=[''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY',k=10)) for _ in range(40)]
            def write(p, seqs):
                p.write_text('sequence\tcharge\tim\n'+''.join(f'{s}\t{z}\t{1+i*.01+z*.15}\n' for i,s in enumerate(seqs) for z in (2,3)))
            write(d/'tr.tsv',seqs[:24]);write(d/'va.tsv',seqs[24:32]);write(d/'test.tsv',seqs[32:])
            cfg={'threads':1,'batch_size':16,'architecture':{'emb':4,'ch':4,'blocks':[1],'head':8}}
            train(d/'tr.tsv',d/'model.uimc',units='inverse_mobility',method='test',from_scratch=True,validation=d/'va.tsv',epochs=2,config=cfg)
            m,info,meta=load_model(d/'model.uimc',target='im')
            with self.assertRaises(ValueError):load_model(d/'model.uimc')
            self.assertEqual(info['n_global'],12)
            self.assertEqual(meta['charge_range'],[2,3])
            # Artifact family and feature width must not be silently reinterpreted.
            import struct
            blob=(d/'model.uimc').read_bytes()
            for bad in (b'URTC'+blob[4:], blob[:4]+struct.pack('<I',2)+blob[8:], blob+b'extra', blob[:60]):
                (d/'corrupt.uimc').write_bytes(bad)
                with self.assertRaises(ValueError): load_model(d/'corrupt.uimc',target='im')
            # Fine-tuning shares the parent-lineage guard and epoch-zero baseline.
            with self.assertRaisesRegex(ValueError,'starting model'):
                train(d/'tr.tsv',d/'bad.uimc',units='inverse_mobility',method='adapt',init=d/'model.uimc',validation=d/'va.tsv',epochs=1,config={'threads':1})
            report=train(d/'tr.tsv',d/'adapt.uimc',units='inverse_mobility',method='adapt',init=d/'model.uimc',validation=d/'test.tsv',epochs=1,config={'threads':1})
            self.assertLessEqual(report['validation']['mae'],report['baseline']['mae']+1e-7)
            exe=os.environ.get('UNO_TEST_BINARY')
            if not exe:return
            def run(*args,ok=True):
                r=subprocess.run([exe,'im',*map(str,args)],capture_output=True,text=True)
                self.assertEqual(r.returncode==0,ok,r.stderr)
                return r
            import json, sys
            (d/'config.json').write_text(json.dumps(cfg))
            run('train','--data',d/'tr.tsv','--out',d/'cli.uimc','--units','inverse-mobility','--method','CLI test','--from-scratch','--epochs',1,'--config',d/'config.json','--python',sys.executable)
            run('inspect',d/'cli.uimc')
            run('predict','--model',d/'model.uimc','--data',d/'test.tsv','--out',d/'pred.tsv')
            rows=list(csv.DictReader((d/'pred.tsv').open(),delimiter='\t'))
            queries=[(r['sequence'],int(r['charge'])) for r in rows]
            oracle=predict(m,queries,info)
            np.testing.assert_allclose([float(r['predicted_im']) for r in rows],oracle,atol=2e-6,rtol=2e-6)
            run('evaluate','--model',d/'model.uimc','--data',d/'test.tsv','--units','inverse-mobility')
            run('evaluate','--model',d/'model.uimc','--data',d/'tr.tsv','--units','inverse-mobility',ok=False)
            run('evaluate','--model',d/'model.uimc','--data',d/'test.tsv','--units','inverse-mobility','--calibration-data',d/'test.tsv',ok=False)
            run('predict','--model',d/'model.uimc','--data',d/'test.tsv','--out',d/'pred.tsv',ok=False)

if __name__=='__main__':unittest.main()
