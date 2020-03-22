#!/usr/bin/env python

import sys
import parmed as pmd
import numpy as np
from scipy.spatial import distance

if len(sys.argv) < 2:
	print "Usage: molecular_diameter.py <mymolecule.mol2>"
	exit(1)

mol  = pmd.load_file(sys.argv[1])
crds = mol.coordinates

dist = distance.cdist(crds, crds, 'euclidean')
print np.max(dist)
exit(0)