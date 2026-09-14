#!/usr/bin/env python3
"""Execution-only concurrency amendment: IonStar alongside the draining Astral queue.
No new settings, cohorts, seeds, or output routes. Run the original driver after
both queues finish to verify reuse and execute the dedicated serialized timings.
"""
from matrix_release_study import *
if __name__=='__main__':
    dataset='ionstar'
    parallel([(extract,(dataset,r)) for r in rows(dataset)])
    parallel([(align,(dataset,n,rr)) for n,rr in cohorts(dataset).items()])
    tasks=[]
    for seed in [0,1,2]:
        tasks += [(search,(dataset,'single_'+r['run'],seed,'single',OUT/dataset/'koth'/r['run']/'features.tsv')) for r in rows(dataset)]
        tasks += [(search,(dataset,arm+'_'+n,seed,arm,OUT/dataset/'aligned'/n)) for n in cohorts(dataset) for arm in ['joint','joint_cell']]
    parallel(tasks)
    write(OUT/'ionstar-queue-complete.json',dict(searches=78))
