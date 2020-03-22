#!/usr/bin/env python

import sys
import numpy as np
from scipy.spatial import distance
import parmed as pmd

if len(sys.argv) != 4:
   print "Usage: filter_clashes.py <target.pdb> <query-filter-structure.pdb> <cutoff-distance> "
   exit(1)

mol_ref_path = sys.argv[1]
mol_cut_path = sys.argv[2]
cutoff  = sys.argv[3]
cutoff  = float(cutoff)**2

mol_ref = pmd.load_file(mol_ref_path)["!@H="]
mol_cut = pmd.load_file(mol_cut_path)["!@H="]

dists_ref  = distance.cdist(mol_ref.coordinates, mol_cut.coordinates, metric="sqeuclidean")
dists_self = distance.cdist(mol_cut.coordinates, mol_cut.coordinates, metric="sqeuclidean")

filter_ref  = np.where(dists_ref<cutoff)
filter_self = np.where(dists_self<cutoff)

filter_list  = list()
filter_block = list()

### Filter everything that crashes with reference
for res in filter_ref[1]:
    if res not in filter_list:
        filter_list.append(res)
        filter_block.append(res)

### Filter everything that crashes with itself
for res0, res1 in zip(filter_self[0], filter_self[1]):
    ### Before checking for some special cases, 
    ### do some general checks.
    if res0 == res1:
        continue
    if res0 in filter_list \
    or res0 in filter_block:
        continue
    if res1 in filter_list \
    or res1 in filter_block:
        continue
    ### Case 1: res0 crashes with reference but res1 not
    if res0 in filter_list \
    and not res1 in filter_list:
        filter_block.append(res1)
    ### Case 2: res1 crashes with reference but res0 not
    if res1 in filter_list \
    and not res0 in filter_list:
        filter_block.append(res0)
    ### Case 3: res0 and res1 crash with reference
    if res0 in filter_list \
    and res1 in filter_list:
        filter_block.append(res0)
        filter_block.append(res1)
    ### Do filtering
    if not res0 in filter_block \
    and not res0 in filter_list:
        filter_list.append(res0)
        filter_block.append(res0)
    if not res1 in filter_block \
    and not res1 in filter_list:
        filter_list.append(res1)
        filter_block.append(res1)

### Strip off all residues that should be filtered
N_filter  = len(filter_list)
if N_filter>0:
    strip_str = ':'
    for i in range(N_filter):
        res = filter_list[i]
        strip_str += '%d' %res
        if i != N_filter-1:
            strip_str += ','
    mol_cut.strip(strip_str)
mol_cut.write_pdb("filter.pdb")