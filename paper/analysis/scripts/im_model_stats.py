"""Derived SI values for the frozen IM training and search experiment."""
import json
from pathlib import Path

BENCH = Path(__file__).resolve().parents[3] / 'benchmark'


def add(st):
    report = json.loads((BENCH / 'im_cnn_20260914.json').read_text())
    diag = report['prediction_diagnostics']
    train = diag['training_report']
    for name, value in [('train_obs', train['train_observations']), ('validation_obs', train['validation_observations']), ('test_obs', diag['public_test']['observations']), ('best_epoch', train['best_epoch']), ('order_n', diag['order_diagnostic']['n'])]:
        st.add('im.'+name, value, fmt=',d', desc='IM '+name.replace('_', ' '), between=(1, 1000000))
    for name, value in report['verification']['sequence_disjoint_splits'].items():
        st.add('im.'+name+'_seq', value, fmt=',d', desc='IM '+name+' disjoint sequences', between=(1, 1000000))
    for name, value in [('mae', diag['public_test']['mae']), ('rmse', diag['public_test']['rmse']), ('order_original', diag['order_diagnostic']['original']['mae']), ('order_reversed', diag['order_diagnostic']['internal_sequence_reversed']['mae'])]:
        st.add('im.'+name, value, fmt='.5f', desc='IM proportional inverse mobility '+name, between=(0, 1))
    st.add('im.r2', diag['public_test']['r2'], fmt='.4f', desc='IM independent public test R2', between=(0, 1))
    for name in ('improved', 'worse', 'tied'):
        st.add('im.'+name, report['study']['outcomes']['cnn_soft'][name], fmt='d', desc='CNN versus off paired searches '+name, between=(0, 9))
