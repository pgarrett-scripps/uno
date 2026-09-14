import tempfile
import csv
from pathlib import Path
import unittest
import numpy as np
from matrix_study import fasta_records
from matrix_report import accession_union, evaluate, list_counts, quant_metrics, compare_cohort
from matrix_evidence import load_evidence


class MatrixStudyTests(unittest.TestCase):
    def test_fasta_description_can_contain_angle_bracket(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'input.fasta'
            p.write_text('>sp|A|A_HUMAN description > something\nABC\nDEF\n>sp|B|B_YEAST\nGHI\n')
            self.assertEqual(list(fasta_records(p)), [
                ('sp|A|A_HUMAN description > something', 'ABCDEF'), ('sp|B|B_YEAST', 'GHI')])

    def test_cohort_union_deduplicates_expanded_groups(self):
        a = {'q': .001, 'species': 'HUMAN', 'members': ['A_HUMAN', 'B_HUMAN']}
        b = {'q': .001, 'species': 'HUMAN', 'members': ['B_HUMAN', 'C_HUMAN']}
        y = {'q': .001, 'species': 'YEAST', 'members': ['Y_YEAST']}
        union = accession_union([[a, y], [b, y]])
        self.assertEqual(len(union), 4)
        self.assertEqual(list_counts(union, 3)['entrapment_fdp'], 1)

    def test_entrapment_ties_are_conservative_and_order_independent(self):
        rows = [dict(accession='H', score=10., species='HUMAN', q=.001),
                dict(accession='Y', score=10., species='YEAST', q=.001)]
        a, b = evaluate(rows, 1), evaluate(rows[::-1], 1)
        self.assertEqual(a, b)
        self.assertEqual(a['entrapment_depth']['0.01'], 0)
        self.assertEqual(a['internal_1pct']['entrapment_fdp'], 1)

    def evidence_fixture(self, directory, mismatch=False):
        meta = ['massCalib','mz','charge','rtApex','im','combined_score','seed_run','n_contributing_runs']
        vals = ['1000','501','2','10','','0.9','A','1']
        for filename, cells in [('intensity_matrix.tsv',['100','200']),('qvalue_matrix.tsv',['0.0050','0.5000'])]:
            with open(directory/filename,'w') as f:
                w=csv.writer(f,delimiter='\t');w.writerow(meta+['A','B']);w.writerow(vals+cells)
        details = ['feature_idx','massCalib','mz','charge','refRtApex','refIm','seed_combined_score',
                   'seed_run','n_contributing_runs','run_name','is_decoy','is_mbr','intensity',
                   'q_value','spectral_bhattacharyya','coelution','n_isotopes_found','ppm_error','rt_diff']
        with open(directory/'lfq_details.tsv','w') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(details)
            w.writerow(['0']+vals+['A','false','false','100','0.005','0.9','0.8','3','-1','-0.2'])
            changed=vals.copy()
            if mismatch: changed[0]='1001'
            w.writerow(['0']+changed+['B','false','true','200','0.5','0.2','0.1','1','2','0.4'])

    def test_evidence_separates_confident_cells_from_mbr(self):
        with tempfile.TemporaryDirectory() as d:
            directory=Path(d);self.evidence_fixture(directory)
            _,_,_,mask,quality,summary,_=load_evidence(directory)
            self.assertEqual(mask.tolist(),[[True,False]])
            self.assertEqual(quality[0,:3].tolist(),[.5,.5,0.])
            self.assertEqual(summary['confident_direct_cells'],1)
            self.assertEqual(summary['confident_mbr_cells'],0)

    def test_evidence_rejects_same_length_but_mismatched_sidecar(self):
        with tempfile.TemporaryDirectory() as d:
            directory=Path(d);self.evidence_fixture(directory,mismatch=True)
            with self.assertRaisesRegex(ValueError,'Metadata mismatch'):
                load_evidence(directory)

    def test_sample_support_requires_confident_same_run_evidence(self):
        from matrix_astral_report import precursor_support
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'lfq.tsv'
            p.write_text('peptide\tcharge\tdecoy\tproteins\tmass\trt\tA\tB\n'
                         'PEPTIDE\t2\t0\tsp|P1|P1_HUMAN\t1000.00000\t10.0000\t100\t200\n'
                         'DECOY\t2\t1\tsp|P1|P1_HUMAN\t1000.00000\t10.0000\t100\t200\n')
            cells={'run_names':np.array(['A','B']),
                   'confident':np.array([[True,False]]),
                   'direct':np.array([[False,True]])}
            cache={'all':(cells,{('1000.00000','10.0000','2'):0})}
            msms={'A':{'precursors':['PEPTIDE/2']},'B':{'precursors':['PEPTIDE/2']}}
            result=precursor_support(p,'all',{'P1'},msms,cache)
            self.assertEqual(result['mapped_assignment_rows'],1)
            self.assertEqual(result['runs']['A']['confident_mbr_rows'],1)
            self.assertEqual(result['runs']['A']['msms_supported_confident_precursors'],1)
            self.assertEqual(result['runs']['B']['positive_assignment_rows'],1)
            self.assertEqual(result['runs']['B']['unique_confident_precursors'],0)

    def test_known_ratio_metrics_do_not_impute_missing_cells(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'abundance.tsv'
            runs=[f'human_ecoli_{c}_rep{n}' for c in 'AE' for n in range(4)]
            with open(p,'w') as f:
                w=csv.writer(f,delimiter='\t')
                w.writerow(['accession','n_group_members','group_members','n_peptides']+runs)
                for n in range(3):w.writerow([f'H{n}_HUMAN',1,f'H{n}_HUMAN',2]+[10]*8)
                w.writerow(['E_ECOLI',1,'E_ECOLI',2]+[10]*4+[30]*4)
                w.writerow(['M_ECOLI',1,'M_ECOLI',2]+[10,10,0,0,30,30,0,0])
            result=quant_metrics(p)
            e=result['species']['ECOLI']
            self.assertAlmostEqual(e['median_abs_log2_ratio_error'],0)
            self.assertEqual(e['median_replicate_cv'],0)
            self.assertEqual(e['positive_cell_fraction'],.75)
            common=quant_metrics(p,fixed_factors=result['size_factors_log2'],common_groups={('M_ECOLI',)})
            self.assertEqual(common['species']['ECOLI']['groups'],1)
            self.assertEqual(common['species']['ECOLI']['positive_cell_fraction'],.5)
            self.assertEqual(common['size_factors_log2'],result['size_factors_log2'])

    def test_added_ids_exclude_entrapment_and_keep_per_run_scope(self):
        def row(sp,members):return {'q':.001,'species':sp,'members':members}
        joint=[row('HUMAN',['H1_HUMAN','H2_HUMAN']),row('YEAST',['Y_YEAST'])]
        singles={'A':[row('HUMAN',['H2_HUMAN']),row('ECOLI',['E_ECOLI'])],'B':[]}
        result=compare_cohort(joint,singles)
        self.assertEqual(result['joint_only_accessions'],['H1_HUMAN'])
        self.assertEqual(result['separate_only_accessions'],['E_ECOLI'])
        self.assertEqual(result['intersection'],1)
        self.assertEqual(result['per_run']['A']['fraction_recovered'],.5)
        self.assertIsNone(result['per_run']['B']['fraction_recovered'])


if __name__ == '__main__': unittest.main()
