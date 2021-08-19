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
parser.add_argument('--lattice_factors_x', '-lx',type=int,\
                                                 help='Repeat unit cell in x direction from a to b',\
                                                 default=[-1,1],\
                                                 nargs=2)
parser.add_argument('--lattice_factors_y', '-ly',type=int,\
                                                 help='Repeat unit cell in y direction from a to b',\
                                                 default=[-1,1],\
                                                 nargs=2)
parser.add_argument('--lattice_factors_z', '-lz',type=int,\
                                                 help='Repeat unit cell in z direction from a to b',\
                                                 default=[-1,1],\
                                                 nargs=2)
parser.add_argument('--output',          '-o',   type=str,\
                                                 help='Output filename',\
                                                 default='lattice.pdb')
parser.add_argument('--verbose',         '-v',   type=int,\
                                                 help='Verbosity mode',\
                                                 default=0)

def symexpcell(
    prefix='mate',
    pymol_object=None, 
    reference=None, 
    a=0, 
    b=0, 
    c=0):

    '''
DESCRIPTION

    Adapted from supercell.py: https://pymolwiki.org/index.php/Supercell

    Creates all symmetry-related objects for the specified object that
    occur with their bounding box center within the unit cell.

USAGE

    symexpcell prefix, object, [a, b, c]

ARGUMENTS

    prefix = string: prefix of new objects

    pymol_object = string: pymol object for which to create symmetry mates

    a, b, c = integer: create neighboring cell {default: 0,0,0}

SEE ALSO

    symexp, http://www.pymolwiki.org/index.php/SuperSym
    '''
    if pymol_object == None:
        pymol_object = pymol.cmd.get_object_list()[0]
    if reference == None:
        reference = pymol_object

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

    matrices    = pymol.xray.sg_sym_to_mat_list(spacegroup)
    N_matrices  = len(matrices)
    crds_list   = np.zeros(
        (N_matrices,
         pymol.cmd.count_atoms(pymol_object),
         3
         )
        )
    for mat_i in range(N_matrices):
        mat   = np.matrix(matrices[mat_i])
        shift = np.floor(mat * center_cell)
        mat[0:3,3] -= shift[0:3,0]
        mat[0:3,3] += extra_shift

        mat = basis * mat * basis.I
        mat_list = list(mat.flat)

        name = '%s%d' % (prefix, mat_i)
        pymol.cmd.create(name, pymol_object)
        pymol.cmd.transform_object(name, mat_list, 0)
        crds_list[mat_i] = pymol.cmd.get_coords(name, 1)

    unique_objects = [0]
    for mat_i in range(N_matrices):
        if mat_i in unique_objects:
            continue
        ### Compute distance between molecule generated through
        ### matrix `mat_i` and all other molecules
        dists = np.linalg.norm(
                crds_list[unique_objects]-crds_list[mat_i],
                axis=2
                )
        too_close_list = np.where(dists < 1e-10)[0]
        if too_close_list.size == 0:
            unique_objects.append(mat_i)

    for obj_i in range(len(unique_objects)):
        unique_objects[obj_i] = '%s%d' % (prefix, unique_objects[obj_i])

    return unique_objects


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
    pymol.cmd.set("retain_order", "1")

    pymol.cmd.load(args.reference, "reference")
    object_count = 0

    for input_idx, input in enumerate(args.input):
        if verbose:
            print("Loading %s..." %input)
        pymol.cmd.load(input, "input")
        for a in range(args.lattice_factors_x[0],args.lattice_factors_x[1]+1):
            for b in range(args.lattice_factors_y[0],args.lattice_factors_y[1]+1):
                for c in range(args.lattice_factors_z[0],args.lattice_factors_z[1]+1):
                    if verbose:
                        print("Generating unit cell",a,b,c)
                    object_list = symexpcell("i_", "input", "reference", a, b, c)
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
