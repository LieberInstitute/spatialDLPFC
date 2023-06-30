#   The o
import scanpy as sc
import session_info
from pyhere import here
from pathlib import Path

orig_ann_path = Path(
    here(
        'raw-data', 'psychENCODE', 'version6', 'SZBDMulti-Seq_annotated.h5ad'
    )
)
out_dir = Path(
    here('raw-data', 'psychENCODE', 'version6', 'SZBDMulti-Seq_subset')
)

out_dir.mkdir(exist_ok = True)

ann = sc.read(orig_ann_path)

#   Explore how much data is from each diagnosis
ann.obs['individualID'].str.fullmatch(r'^BD.*').sum()
ann.obs['individualID'].str.fullmatch(r'^CON.*').sum()
ann.obs['individualID'].str.fullmatch(r'^SZ.*').sum()

for diagnosis in ('CON', 'BD', 'SZ'):
    #   Subset to this diagnosis and write
    this_ann = ann[ann.obs['individualID'].str.fullmatch(rf'^{diagnosis}.*'), :]
    this_ann.write(out_dir / f'SZBDMulti-Seq_annotated_{diagnosis}.h5ad')

session_info.show()
