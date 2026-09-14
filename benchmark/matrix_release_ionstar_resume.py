#!/usr/bin/env python3
"""Resume the identical IonStar search inventory with four workers after Astral.
Completed job receipts are verified and reused. Interrupted attempts are archived.
"""
from concurrent.futures import ThreadPoolExecutor,as_completed
from matrix_release_study import OUT,rows,cohorts,search,write
if __name__=='__main__':
    dataset='ionstar';tasks=[]
    for seed in [0,1,2]:
        tasks += [(dataset,'single_'+r['run'],seed,'single',OUT/dataset/'koth'/r['run']/'features.tsv') for r in rows(dataset)]
        tasks += [(dataset,arm+'_'+n,seed,arm,OUT/dataset/'aligned'/n) for n in cohorts(dataset) for arm in ['joint','joint_cell']]
    assert len(tasks)==78
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(search,*args) for args in tasks]
        for future in as_completed(futures):future.result()
    write(OUT/'ionstar-queue-complete.json',dict(searches=78))
