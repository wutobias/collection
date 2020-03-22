#!/usr/bin/env python
"""

Calculate center of geometry for given residue(s)
in PDB file.
Usage: get_cog.py <my_protein.pdb> #RESID #RESID #RESID #RESID ...

"""

import os
import sys
import numpy as np
import parmed as pmd

if len(sys.argv) < 3:
	print "Usage: get_cog.py <my_protein.pdb> #RESID #RESID #RESID #RESID ..."
	exit(0)

res_list = sys.argv[2:]
#Generate selection
sele = ":"
for r in res_list:

	sele += r
	if r != res_list[-1]:
		sele += ","

pdb      = pmd.load_file(sys.argv[1])

cog = np.mean(pdb[sele].coordinates, axis=0)

print "%6.3f" %(cog[0]), "%6.3f" %(cog[1]), "%6.3f" %(cog[2])