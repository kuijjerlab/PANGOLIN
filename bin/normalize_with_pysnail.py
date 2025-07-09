# normalize_with_pysnail.py
import os
from pysnail import Dataset, qsmooth

xprs = os.path.realpath(snakemake.input["counts"])
groups = os.path.realpath(snakemake.input["groups"])
output_file = os.path.realpath(snakemake.output[0])

dataset = Dataset(xprs, groups, **{'index_col': 0, 'sep': '\t'})
xprs_norm, qstat = qsmooth(dataset, aggregation='auto', threshold=0.2)
xprs_norm.to_csv(output_file, sep='\t')