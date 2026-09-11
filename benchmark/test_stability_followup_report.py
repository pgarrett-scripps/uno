"""Scientific gate fixtures plus synthetic-executable provenance integration."""
from pathlib import Path
import contextlib,copy,importlib.util,io,json,tempfile,unittest
from unittest.mock import patch
import stability_followup_report as report


def fixture():
 rows=[]
 for ds in sorted(report.DATASETS['lowinput']):
  for seed in report.SEEDS['screening']['lowinput']:
   for variant in ('baseline','candidate'):
    rows.append(dict(dataset=ds,seed=seed,variant=variant,depth_1pct=100,depth_5pct=100,rows=200,internal_q1_targets=100,internal_q1_human=100,internal_q1_yeast=0,internal_q1_other=0,fdp=0.,entrapment1_human_set={f'H{i}' for i in range(100)},reported_q1_target_set={f'H{i}' for i in range(100)}))
 return rows

def summarize(rows):return report.summarize_panel(rows,1.,'lowinput','lowinput-screen')
def checks(result):return {r['criterion']:r['passed'] for r in result['checks']}
def change(records,variant,seed,ds=None,**values):
 ds=ds or sorted(report.DATASETS['lowinput'])[0]
 r=next(r for r in records if (r['dataset'],r['seed'],r['variant'])==(ds,seed,variant));r.update(values);return r

def rescue(rows):
 change(rows,'baseline',0,depth_1pct=0,depth_5pct=0,rows=0,internal_q1_targets=0,internal_q1_human=0,fdp=None,entrapment1_human_set=set(),reported_q1_target_set=set())

class Gates(unittest.TestCase):
 def test_failure_rescue_and_exact_core_stress(self):
  rows=fixture();rescue(rows);r=summarize(rows)
  self.assertTrue(all(checks(r).values()));self.assertEqual(r['rescued_only']['all']['n'],1)
  self.assertEqual(r['aggregates']['core']['baseline']['n'],18);self.assertEqual(r['aggregates']['stress']['baseline']['n'],12)
  self.assertEqual(r['rescued_only']['core']['n'],1);self.assertEqual(r['rescued_only']['stress']['n'],0)
 def test_new_empty_rejected_even_with_rescue(self):
  rows=fixture();rescue(rows);change(rows,'candidate',1,internal_q1_targets=0)
  self.assertFalse(checks(summarize(rows))['lowinput_no_new_empty_q1'])
 def test_per_file_median_not_masked_by_pooled_successes(self):
  rows=fixture();rescue(rows)
  ds=sorted(report.DATASETS['lowinput'])[-1]
  for seed in report.SEEDS['screening']['lowinput']:change(rows,'candidate',seed,ds,depth_1pct=94)
  r=summarize(rows);self.assertFalse(checks(r)['lowinput_no_file_median_loss_over5pct'])
  self.assertEqual(r['aggregates']['all']['candidate']['median_file_depth_1pct'],100)
  for seed in report.SEEDS['screening']['lowinput']:change(rows,'candidate',seed,ds,depth_1pct=95)
  self.assertTrue(checks(summarize(rows))['lowinput_no_file_median_loss_over5pct'])
 def test_fdp_boundary_and_above(self):
  rows=fixture()
  for r in rows:r.update(internal_q1_human=200,internal_q1_targets=200+(r['variant']=='candidate'),internal_q1_yeast=int(r['variant']=='candidate'),fdp=.005 if r['variant']=='candidate' else 0.)
  self.assertTrue(checks(summarize(rows))['lowinput_pooled_fdp_delta_at_most0.5pp'])
  for r in rows:
   if r['variant']=='candidate':r.update(internal_q1_yeast=2,internal_q1_targets=202,fdp=.01)
  self.assertFalse(checks(summarize(rows))['lowinput_pooled_fdp_delta_at_most0.5pp'])
 def test_single_over10_loss_is_disclosure_only(self):
  rows=fixture();rescue(rows);change(rows,'candidate',1,depth_1pct=89)
  r=summarize(rows);self.assertTrue(all(checks(r).values()));self.assertEqual(len(r['losses_over10pct']),1)
 def test_variability_improvement_without_rescue(self):
  rows=fixture()
  for r in rows:
   if r['variant']=='baseline':r['depth_1pct']={0:50,1:100,2:150,5:100,8:100}[r['seed']]
  r=summarize(rows);self.assertTrue(checks(r)['lowinput_improves_failures_or_variability']);self.assertEqual(r['relative_ranges']['all']['candidate'],0.)
  self.assertFalse(checks(summarize(fixture()))['lowinput_improves_failures_or_variability'])
 def test_undefined_fdp_ranges_and_empty_empty_overlaps(self):
  rows=fixture()
  for r in rows:r.update(depth_1pct=0,depth_5pct=0,rows=0,internal_q1_human=0,internal_q1_targets=0,fdp=None,entrapment1_human_set=set(),reported_q1_target_set=set())
  r=summarize(rows);self.assertIsNone(checks(r)['lowinput_pooled_fdp_delta_at_most0.5pp']);self.assertIsNone(r['relative_ranges']['all']['baseline']);self.assertEqual(len(r['relative_ranges']['all']['excluded']['baseline']),6)
  self.assertTrue(all(p['paired_reported_q1_target_jaccard'] is None and p['paired_entrapment1_human_jaccard'] is None for p in r['pairs']))
  self.assertTrue(all(p['reported_q1_target_jaccard'] is None and p['entrapment1_human_jaccard'] is None for p in r['within_variant_overlap']))
 def test_other_and_overlap_semantics(self):
  rows=fixture()
  a=change(rows,'baseline',0,internal_q1_other=2,internal_q1_targets=102,reported_q1_target_set={'OTHER1','OTHER2'},entrapment1_human_set={'H1'})
  b=change(rows,'candidate',0,internal_q1_other=1,internal_q1_targets=101,reported_q1_target_set={'OTHER1'},entrapment1_human_set={'H1'})
  r=summarize(rows);p=next(p for p in r['pairs'] if p['dataset']==a['dataset'] and p['seed']==0)
  self.assertEqual(p['paired_reported_q1_target_jaccard'],.5);self.assertEqual(p['paired_entrapment1_human_jaccard'],1.)
  self.assertEqual(r['aggregates']['all']['baseline']['other'],2);self.assertEqual(r['aggregates']['all']['candidate']['other'],1)
 def test_descriptive_reported_range_and_cv_not_gates(self):
  rows=fixture();rescue(rows);r=summarize(rows);d=r['descriptive_reported_q1_variability']['internal_q1_targets'];self.assertEqual(len(d['common_range_files']),6)
  first=d['files'][0];self.assertEqual(first['baseline']['relative_range'],1.);self.assertGreater(first['baseline']['sample_cv'],0.);self.assertEqual(first['candidate']['sample_cv'],0.)
  for row in rows:row['internal_q1_targets']=0;row['internal_q1_human']=0
  d=report.reported_variability(rows)['internal_q1_targets'];self.assertEqual(d['common_range_files'],[]);self.assertIsNone(d['median_relative_range']['baseline']);self.assertIsNone(d['files'][0]['baseline']['sample_cv'])
 def test_reserved_seed_membership(self):
  self.assertEqual(report.SEEDS['confirmation'],{'lowinput':[9,10],'primary':[11,12],'raw':[5,6]})

class Provenance(unittest.TestCase):
 def test_synthetic_search_metadata_and_rejections(self):
  # These executables emit three fixed protein rows; no real UNO runs.
  with tempfile.TemporaryDirectory(prefix='followup-report-fixture-',dir='/tmp') as tmp:
   root=Path(tmp);helper=report.helper
   def stub(name):
    p=root/name;p.write_text('#!/usr/bin/env python3\nfrom pathlib import Path\nimport sys\nseed=sys.argv[sys.argv.index("--seed")+1]\np=Path(sys.argv[sys.argv.index("--out-proteins")+1])\np.write_text("accession\\tdecoy\\tbinomial_score\\tq_value\\nA_HUMAN\\t0\\t10\\t0.02\\nB_HUMAN\\t0\\t9\\t0.005\\n'+name+'"+seed+"_MOUSE\\t0\\t8\\t0.005\\n")\nprint("synthetic fixture only")\n');p.chmod(0o755);return p
   baseline=stub('baseline.py');candidate=stub('candidate.py');features=root/'features';features.write_text('synthetic features')
   model=root/'model';model.write_text('synthetic model');config=root/'config';config.write_text('synthetic config');inventory=root/'inventory';inventory.write_text('synthetic inventory')
   fasta=root/'fasta';fasta.write_text(''.join(f'>H{i}_HUMAN\nAA\n' for i in range(20647))+''.join(f'>Y{i}_YEAST\nAA\n' for i in range(6066)))
   specs={ds:dict(dataset=ds,panel='lowinput',astral_stem='',features=str(features)) for ds in report.DATASETS['lowinput']};out=root/'panel'
   with patch.object(helper,'inventory',return_value=(specs,{})),patch.object(helper,'INVENTORY',inventory),contextlib.redirect_stdout(io.StringIO()):
    rc=helper.main(['--baseline',str(baseline),'--candidate',str(candidate),'--output',str(out),'--datasets','lowinput','--seeds','0,1,2,5,8','--held-out-seeds','9-10','--threads','4','--jobs','2','--env','UNO_RT_CACHE_DIR=/tmp/uno-stability/rt-cache','--fasta',str(fasta),'--model',str(model),'--config',str(config)])
   self.assertEqual(rc,0)
   with self.assertRaisesRegex(ValueError,'wrong engine'):report.load_panel(out,'lowinput',[0,1,2,5,8],report.sha(candidate),{})
   with patch.object(report,'V3_SHA',report.sha(baseline)):
    rows,scale,tables,_=report.load_panel(out,'lowinput',[0,1,2,5,8],report.sha(candidate),{})
    self.assertEqual(len(rows),60);self.assertTrue(all(r['internal_q1_other']==1 for r in rows))
    r=report.summarize_panel(rows,scale,'lowinput','lowinput-screen',tables)
    self.assertTrue(all(p['paired_entrapment1_human_jaccard']==1 and p['paired_reported_q1_target_jaccard']==1/3 for p in r['pairs']))
    p=out/'runs'/f"{sorted(report.DATASETS['lowinput'])[0]}.seed0.baseline.tsv";saved=p.read_bytes();p.write_bytes(saved+b'tampered\n')
    with self.assertRaisesRegex(ValueError,'hash mismatch'):report.load_panel(out,'lowinput',[0,1,2,5,8],report.sha(candidate),{})
    p.write_bytes(saved)
    # A newly certified but incomplete summary must still fail exact coverage.
    summary=out/'runs.tsv';original=summary.read_bytes();summary.write_text('\n'.join(summary.read_text().splitlines()[:-1])+'\n')
    cert_path=out/'summary_provenance.json';saved_cert=cert_path.read_bytes();cert=json.loads(saved_cert);cert['summaries']['runs.tsv']=report.sha(summary);cert_path.write_text(json.dumps(cert))
    with self.assertRaisesRegex(ValueError,'coverage'):report.load_panel(out,'lowinput',[0,1,2,5,8],report.sha(candidate),{})
    summary.write_bytes(original);cert_path.write_bytes(saved_cert)
    with contextlib.redirect_stdout(io.StringIO()):
     report.main(['--stage','lowinput-screen','--arm','folds','--candidate-sha256',report.sha(candidate),'--lowinput',str(out),'--output',str(root/'result')])
    result=json.loads((root/'result.json').read_text());self.assertFalse(result['gate_pass']);self.assertFalse(result['eligible_for_promotion'])
    self.assertEqual(result['panels']['lowinput']['aggregates']['all']['baseline']['other'],30)
    # A rejected prior screen cannot be rescued by unseen-only evidence.
    prior=root/'rejected-screen.json';failed=copy.deepcopy(result);failed['stage']='screening';prior.write_text(json.dumps(failed))
    with self.assertRaisesRegex(ValueError,'passing matching frozen screening'):
     report.main(['--stage','confirmation','--arm','folds','--candidate-sha256',report.sha(candidate),'--lowinput',str(out),'--primary',str(out),'--raw',str(out),'--screening-report',str(prior),'--output',str(root/'forbidden')])
    # No passing screening evidence means no unseen-based promotion.
    with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
     report.main(['--stage','confirmation','--arm','folds','--candidate-sha256',report.sha(candidate),'--lowinput',str(out),'--primary',str(out),'--raw',str(out),'--output',str(root/'bad')])

if __name__=='__main__':unittest.main(verbosity=2)
