#!/usr/bin/env python

import numpy as np
from scipy import spatial
import parmed as pmd
import sys


if len(sys.argv) != 3:

    print "Usage: normalize_B.py <my_protein.PDB> <radius> [Angstrom]"
    exit(0)

if sys.argv[1].endswith("pdb"):

    name     = sys.argv[1].replace(".pdb", "")

elif sys.argv[1].endswith("PDB"):

    name     = sys.argv[1].replace(".PDB", "")

else:

    name     = sys.argv[1]

pdb      = pmd.load_file(sys.argv[1])

norm_pdb = pmd.load_file(sys.argv[1])

std_pdb  = pmd.load_file(sys.argv[1])

r        = float(sys.argv[2])

print "Input:", sys.argv[1]
print "radius", r, "Ang"

matrix_idx = np.where( spatial.distance.cdist(pdb.coordinates, pdb.coordinates, metric='euclidean') < r )

for idx_1 in range(len(pdb.atoms)):

    if pdb[idx_1].bfactor == 0.0:

        continue
        
    new_b = list()

    neighbor_list = pdb[matrix_idx[1][np.where(matrix_idx[0] == idx_1)]].atoms

    for atom_2 in neighbor_list:

        if atom_2.bfactor == 0.0:

            continue

        atom2_name = atom_2.name

        if (atom2_name.startswith("H") \
        or atom2_name.startswith("h")) \
        and atom2_name != "HG":
    
            continue
            
        new_b.append(atom_2.occupancy+atom_2.bfactor)

    norm_pdb[idx_1].bfactor /= np.array(new_b).mean()

    std_pdb[idx_1].bfactor = (std_pdb[idx_1].bfactor - np.array(new_b).mean()) / np.array(new_b).std()

pmd.write_PDB(norm_pdb, name + "_norm.pdb")

pmd.write_PDB(std_pdb, name + "_std.pdb")