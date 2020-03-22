#!/usr/bin/env python

import parmed as pmd
import random
import sys
import os

if not (len(sys.argv) == 5 or len(sys.argv) == 6):

	print "Usage: strip_random.py [-O] <prmtop> <inpcrd> <first_strippable_water> <remaining_wats>"
	exit(1)

overwrite = False
arg_shift = 0
if sys.argv[1] == "-O":
	overwrite = True
	arg_shift = 1
prmtop         = sys.argv[1+arg_shift]
inpcrd         = sys.argv[2+arg_shift]
first_wat      = int(sys.argv[3+arg_shift])
remaining_wats = int(sys.argv[4+arg_shift])

p = pmd.load_file(prmtop)
p.load_rst7(inpcrd)
wats = p[":WAT"]

prmtop_watcnt = len(wats.residues)

if prmtop_watcnt <= remaining_wats:

	print "Warning: Number of Remaining wates lower than present waters."
	exit(1)

wat_list = list()
no_strip = list()
# Get water resids
for wat_idx in range(prmtop_watcnt):

	#Add 1 because parmed internal resiude numbering starts at 1
	if wats.residues[wat_idx].number+1 >= first_wat:
		wat_list.append(wats.residues[wat_idx].number+1)

	else:
		no_strip.append(wats.residues[wat_idx].number+1)

random.shuffle(wat_list)
random.shuffle(wat_list)

strip_mask = ":"
rm_mask    = ""
for strip_count, wat_resid in enumerate(wat_list):

	if strip_count < prmtop_watcnt-remaining_wats-1:

		strip_mask += "%s," %wat_resid
		#Assuming that the final unit is called 'complex'
		rm_mask    += "remove complex complex.%s \n" %wat_resid

	else:

		strip_mask += "%s" %wat_resid
		rm_mask    += "remove complex complex.%s \n" %wat_resid

		break

o = open("strip.cpptraj", "w")
o.write("parm %s \n" %prmtop)
o.write("trajin %s \n" %inpcrd)
o.write("strip %s \n" %strip_mask)
if overwrite:
	o.write("trajout %s restart \n" %inpcrd)
	o.write("trajout %s pdb \n"     %inpcrd.replace("inpcrd", "pdb").replace("incrd", "pdb").replace("rst", "pdb").replace("rst7", "pdb"))
else:
	o.write("trajout strip.inpcrd restart \n")
	o.write("trajout strip.pdb pdb \n")
o.write("go \n")
o.write("clear all \n")
o.write("parm %s \n" %prmtop)
o.write("parmstrip %s \n" %strip_mask)
if overwrite:
	o.write("parmwrite out %s \n" %prmtop)
else:
	o.write("parmwrite out strip.prmtop \n")
o.write("go \n")
o.close()	

if os.path.exists("strip.cpptraj"):

	print "Wrote strip.cpptraj. Please run >> cpptraj -i strip.cpptraj << to complete random water strip."

else:

	print "strip.cpptraj not written. Something went wrong."
	exit(1)

o = open("strip.leap", "w")
o.write(rm_mask)

if not overwrite:
	o.write("saveAmberParm complex strip.prmtop strip.inpcrd \n")
	o.write("savePdb complex strip.pdb \n")
o.close()	

if os.path.exists("strip.leap"):

	print "Wrote strip.leap. Please source strip.lib into your leap file right before you save the prmtop inpcrd files"
	print "Note: cpptraj and leap actions should both have the same effect. Use just one of them, NOT both."
	exit(0)

else:

	print "strip.leap not written. Something went wrong."
	exit(1)
