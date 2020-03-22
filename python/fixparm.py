#!/usr/bin/env python

"""
Sometimes Amber parameter topology files are
corrupt, e.g. containing chain discontinuities.
This python script calls parmed.py and aims to
repair paramter and coordinate file accordingly.
"""

import sys
import os
import parmed as pmd


### Command line checking and preparation
### =====================================

if len(sys.argv) < 3:

	print "Usage: fixparm.py [-O] <mytopology.prmtop> <myrestart.inpcrd>"
	exit(0)

if len(sys.argv) == 3 and sys.argv[1] == "-O":

	print "Usage: fixparm.py [-O] <mytopology.prmtop> <myrestart.inpcrd>"
	exit(0)	

if sys.argv[1] == "-O":

	print "I am in overwrite mode."
	parm = sys.argv[2]
	crd  = sys.argv[3]
	outname = parm.replace(".prmtop", "")

else:

	parm = sys.argv[1]
	crd  = sys.argv[2]
	outname = parm.replace(".prmtop", "")+"_fixed"



### Parse everything to parmed.py using its python API
### ==================================================

### Loading prmtop and inpcrd file into parmed.py..
prmtop = pmd.load_file(parm)
prmtop.load_rst7(crd)

### Repairing...
solvent = pmd.tools.actions.defineSolvent(prmtop,":WAT")
solvent.execute()
Mol     = pmd.tools.actions.defineSolvent(prmtop, True)
Mol.execute()

### Rechecking...
checkit = pmd.tools.actions.checkValidity(prmtop)
checkit.execute()

#### Write out files...
prmtop.write_parm(outname+".prmtop")
prmtop.write_rst7(outname+".inpcrd")
prmtop.write_pdb(outname+".pdb")
