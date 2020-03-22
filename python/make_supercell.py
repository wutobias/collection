#!/usr/bin/env python

import argparse
import numpy as np
import string
import __main__
import xtalmd
from xtalmd.utils import cellbasis

parser = argparse.ArgumentParser(description="Generate crystal lattices by using some reference crystal lattice.\
											  Written by Tobias Wulsdorf, AG Klebe, tobias.wulsdorf@uni-marburg.de")

parser.add_argument('--input',           '-i',   type=str,\
												 help='Input PDB file',\
												 required=True,\
												 nargs='+')
parser.add_argument('--reference',       '-r',   type=str,\
												 help='Reference PDB file. Used for construction of crystal lattice',\
												 required=True)
parser.add_argument('--lattice_factors', '-l',   type=int,\
												 help='Repeat unit cell xyz directions [a,b,c] times.',\
												 default=[1,1,1],\
												 nargs=3)
parser.add_argument('--output',          '-o',   type=str,\
												 help='Output filename',\
												 default='lattice.pdb')
parser.add_argument('--verbose',         '-v',   type=int,\
												 help='Verbosity mode',\
												 default=0)

def symexpcell(prefix='mate', object=None, reference=None, a=0, b=0, c=0):
	'''
DESCRIPTION

	Adapted from supercell.py: https://pymolwiki.org/index.php/Supercell

    Creates all symmetry-related objects for the specified object that
    occur with their bounding box center within the unit cell.

USAGE

    symexpcell prefix, object, [a, b, c]

ARGUMENTS

    prefix = string: prefix of new objects

    object = string: object for which to create symmetry mates

    a, b, c = integer: create neighboring cell {default: 0,0,0}

SEE ALSO

    symexp, http://www.pymolwiki.org/index.php/SuperSym
	'''
	if object == None:
		object = pymol.cmd.get_object_list()[0]
	if reference == None:
		reference = object

	sym = pymol.cmd.get_symmetry(reference)
	cell_edges  = sym[0:3]
	cell_angles = sym[3:6]
	spacegroup  = sym[6]

	basis = cellbasis(cell_angles, cell_edges)
	basis = np.matrix(basis)

	extent = pymol.cmd.get_extent(reference)
	center = sum(np.array(extent)) * 0.5
	center = np.matrix(center.tolist() + [1.0]).T
	center_cell = basis.I * center

	extra_shift = [[float(i)] for i in (a,b,c)]

	i = 0
	matrices    = pymol.xray.sg_sym_to_mat_list(spacegroup)
	object_list = []
	for mat in matrices:
		i += 1

		mat = np.matrix(mat)
		shift = np.floor(mat * center_cell)
		mat[0:3,3] -= shift[0:3,0]
		mat[0:3,3] += extra_shift

		mat = basis * mat * basis.I
		mat_list = list(mat.flat)

		name = '%s%d' % (prefix, i)
		pymol.cmd.create(name, object)
		pymol.cmd.transform_object(name, mat_list, 0)
		object_list.append(name)

	return object_list


if __name__ == "__main__":

	args = parser.parse_args()

	verbose = False
	if args.verbose > 0:
		verbose = True
		
	if verbose:
		__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
	else:
		__main__.pymol_argv = ['pymol','-qQc'] # Pymol: quiet and no GUI

	import pymol
	pymol.finish_launching()

	pymol.cmd.set("pdb_conect_all", "on")
	pymol.cmd.set("connect_mode", "1")

	pymol.cmd.load(args.reference, "reference")
	object_count = 0

	for input_idx, input in enumerate(args.input):
		asymm_count = 1
		if verbose:
			print "Loading %s..." %input
		pymol.cmd.load(input, "input")
		for a in range(-args.lattice_factors[0],args.lattice_factors[0]+1):
			for b in range(-args.lattice_factors[1],args.lattice_factors[1]+1):
				for c in range(-args.lattice_factors[2], args.lattice_factors[2]+1):
#		for a in range(args.lattice_factors[0]+1):
#			for b in range(args.lattice_factors[1]+1):
#				for c in range(args.lattice_factors[2]+1):

					if verbose:
						print "Generating unit cell",a,b,c
					object_list = symexpcell("i_", "input", "reference", a, b, c)
					asymm_count = len(object_list)
					for i in range(len(object_list)):
						chain = string.ascii_uppercase[i]
						pymol.cmd.alter(object_list[i], 'chain="%s"' %chain)
						if verbose:
							pymol.cmd.save("mol%d_%d%d%d_%d_" %(input_idx,a,b,c,i) + args.output, object_list[i])
						pymol.cmd.copy("mol%d" %object_count, object_list[i])
						pymol.cmd.group("global", "mol%d" %object_count)
						object_count += 1
		pymol.cmd.delete("input")

	pymol.cmd.save(args.output, "global")
