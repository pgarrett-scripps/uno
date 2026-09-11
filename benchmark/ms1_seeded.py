#!/usr/bin/env python3
"""Run ms1searchpy with its unseeded randomness locked, to test whether it too
becomes run-to-run reproducible (the fair counterpart of uno's fixed seed).

ms1searchpy already seeds most of its stochastic stages (SEED = 50 for RANSAC,
CV permutation, np.random.seed). Three sources are left unseeded and this
wrapper pins them without editing ms1searchpy:
  1. randomized string hashing, which reorders every `list(set(peptides))`
     before a seeded shuffle consumes it                                   -> PYTHONHASHSEED=0 (re-exec)
  2. the decoy FASTA shuffle (utils.py: module-level `random.shuffle`)      -> random.seed(50)
  3. the LightGBM hyperparameter random search (np.random.RandomState())    -> counter-seeded RandomState
  4. LightGBM training itself (multithreaded, no `seed`/`deterministic`)    -> seed=50, deterministic=True,
                                                                              force_row_wise=True
Everything else (arguments, DeepLC, output) is ms1searchpy as shipped.

  ms1_seeded.py <ms1searchpy args...>
"""
import os, random, sys

# PYTHONHASHSEED must be set BEFORE the interpreter starts, so re-exec once if
# the caller did not set it. This one is load-bearing: ms1searchpy repeatedly
# does `list(set(...))` over peptide strings and then applies a *seeded* shuffle
# to that list (main.py ~1978), so with randomized string hashing the seed
# shuffles a differently-ordered list every run and the whole downstream search
# moves. Seeding `random` alone does not fix it.
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

os.environ.setdefault("TF_DETERMINISTIC_OPS", "1")
random.seed(50)

import numpy as np
_RS = np.random.RandomState
_counter = [50]
class _SeededRandomState(_RS):
    def __init__(self, seed=None):
        if seed is None:
            _counter[0] += 1
            seed = _counter[0]
        super().__init__(seed)
np.random.RandomState = _SeededRandomState

import lightgbm as lgb
_train = lgb.train
def _seeded_train(params, *a, **k):
    params = dict(params)
    params.update(seed=50, deterministic=True, force_row_wise=True)
    return _train(params, *a, **k)
lgb.train = _seeded_train

from ms1searchpy.search import run
sys.exit(run())
