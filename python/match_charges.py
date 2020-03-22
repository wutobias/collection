#!/usr/bin/env python

###############################################################################################
# General Purpose:                                                                            #
#                                                                                             #
# Two molecules A and B, with identical topology and elements,                                #
# but differing coordinates, charges and atom names are                                       #
# matched in order to obtain a molecular representation                                       #
# of a molecule containing cartesian coordinates of molecule A                                #
# and charges and atom names of molecule B.                                                   #
#                                                                                             #
# Usage:                                                                                      #
#                                                                                             #
# match_charges.py <target_molecule.mol2> <template_molecule.mol2> > target_matched.mol2      #
# ...will produce target_matched.mol2 in mol2 format containing                               #
# cartesian coordinates from molecule <target_molecule.mol2>                                  #
# and charges and atom names from <template_molecule.mol2>.                                   #
#                                                                                             #
# T.Wulsdorf, AG Klebe, Marburg. 10/2015                                                      #
#                                                                                             #
###############################################################################################

# Command line checking
# ---------------------
import sys

if len(sys.argv) < 3:
	print "Usage: match_charges.py <target_molecule.mol2> <template_molecule.mol2>"
	exit(0)

# Functions
# ---------

def array2string(array, sep):
   s = ""
   for item in array:
      s += item + sep
   return s


def reorder_it(template, target, match_list):

	file_template = open(template,"r")
	template      = file_template.readlines()
	file_template.close()

	file_target = open(target,"r")
	target      = file_target.readlines()
	file_target.close()



	### read in target molecule coordinates ###
	### =================================== ###

	for index, line in enumerate(target):

		if line == "@<TRIPOS>ATOM\n":
			start_ATOM = index+1

		if line == "@<TRIPOS>BOND\n":
			end_ATOM = index
			break

	target_dict = {}

	for i, line in enumerate(target[start_ATOM:end_ATOM]):

		c              = line.rstrip().split()
		target_dict[i] = [ c[2], c[3], c[4] ]


	### read in template mol2 file ###
	### ========================== ###
	for index, line in enumerate(template):

		if line == "@<TRIPOS>ATOM\n":
			start_ATOM = index+1

		if line == "@<TRIPOS>BOND\n":
			end_ATOM = index
			break

	for i, line in enumerate(template[start_ATOM:end_ATOM]):

		c = line.rstrip().split()
		# Write charges from template to charge culumn c[8] in target mol2
		#
		# template{atomID_in_template:[ atom name, atom charge]}
		c[2] = target_dict[match_list[i]][0]
		c[3] = target_dict[match_list[i]][1]
		c[4] = target_dict[match_list[i]][2]
		template[i+start_ATOM] = "\t" + array2string(c, "\t") + "\n"

	return array2string(template, "")


##################################  M  A  I  N  ##################################################

from rdkit import Chem
from rdkit.Chem import rdFMCS
import random
import os
antechamber="/nfs/progs64/amber/amber14/bin/antechamber"
#antechamber="/usr/local/programs/amber14gnu/bin/antechamber"

target_path   = sys.argv[1]
template_path =	sys.argv[2]

target_tmp   = "/tmp/tmp%d.mol2" %random.randint(0,99999)
template_tmp = "/tmp/tmp%d.mol2" %random.randint(0,99999)

# First we need to make really sure, 
# that both files do have sybyl atomtypes.
# ----------------------------------------

os.system(antechamber + \
		  " -i %s" %target_path + \
		  " -fi mol2" + \
		  " -o %s" %target_tmp + \
		  " -fo mol2" + \
		  " -pf y" + \
		  " -at sybyl > /dev/null" )

os.system(antechamber + \
		  " -i %s" %template_path + \
		  " -fi mol2" + \
		  " -o %s" %template_tmp + \
		  " -fo mol2" + \
		  " -pf y" + \
		  " -at sybyl > /dev/null" )

target        = Chem.MolFromMol2File(target_tmp,   removeHs=False)
template      = Chem.MolFromMol2File(template_tmp, removeHs=False)

# match_list: atom indices in target ordered as template's atoms
match_list    = target.GetSubstructMatches(template)
print reorder_it(template_path, target_path, match_list[0])