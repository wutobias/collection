#!/usr/bin/env python

import argparse
import numpy as np
from scipy.spatial import distance
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
parser.add_argument('--duplicate_cutoff',  '-dc',type=float,\
                                                 help='If any atom between any two molecules are closer than this, one of them molecules be removed.',\
                                                 default=1e-5,\
                                                 nargs=1)
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
    duplicate_cutoff=1e-5,
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

    matrices     = pymol.xray.sg_sym_to_mat_list(spacegroup)
    global_names = pymol.cmd.get_object_list("global")
    N_global     = len(global_names)
    N_matrices   = len(matrices)
    crds_list_mates = np.zeros(
        (N_matrices,
         pymol.cmd.count_atoms(pymol_object),
         3
         )
        )
    crds_list_global = np.zeros(
        (N_global,
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
        crds_list_mates[mat_i] = pymol.cmd.get_coords(name, 1)

    for global_i in range(N_global):
        crds_list_global[global_i] = pymol.cmd.get_coords(global_names[global_i], 1)

    excluded_objects = []
    for mat_i in range(N_matrices):
        if mat_i in excluded_objects:
            continue
        ### Compute distance between molecule generated through
        ### matrix `mat_i` and all other molecules
        for j in range(N_matrices):
            if mat_i == j:
                continue
            dists_mates = distance.cdist(
                    crds_list_mates[j],
                    crds_list_mates[mat_i],
                    )
            too_close_list_mates  = np.where(dists_mates < duplicate_cutoff)[0]
            if too_close_list_mates.size > 0:
                excluded_objects.append(j)
        for j in range(N_global):
            dists_global = distance.cdist(
                    crds_list_global[j],
                    crds_list_mates[mat_i],
                    )
            too_close_list_global = np.where(dists_global < duplicate_cutoff)[0]
            if too_close_list_global.size > 0:
                excluded_objects.append(mat_i)
                break

    unique_objects = list()
    for mat_i in range(N_matrices):
        if not mat_i in excluded_objects:
            unique_objects.append('%s%d' % (prefix, mat_i))
        else:
            pymol.cmd.delete('%s%d' % (prefix, mat_i))

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
    pymol.cmd.set('pdb_use_ter_records', 1)
    pymol.cmd.group("global")

    pymol.cmd.load(args.reference, "reference")
    object_count = 0
    uc_count     = 0

    for input_idx, input in enumerate(args.input):
        if verbose:
            print("Loading %s..." %input)
        pymol.cmd.load(input, "input")
        ### Build the P1 cell
        p1_object_list = symexpcell(
            "p1_cell",
            "input",
            "reference",
            args.duplicate_cutoff,
            0,
            0,
            0)
        
        crds_list = list()
        N_p1_obj  = len(p1_object_list)
        remove_atoms = list()
        for object_name in p1_object_list:
            crds_list.append(pymol.cmd.get_coords(object_name, 1))
        crds_list = np.array(crds_list)
        for object_idx_1 in range(N_p1_obj):
            for object_idx_2 in range(N_p1_obj):
                p1_dists = distance.cdist(
                    crds_list[object_idx_1],
                    crds_list[object_idx_2]
                    )
                valids_1, valids_2 = np.where(p1_dists < args.duplicate_cutoff)
                for i in range(valids_1.size):
                    if valids_1[i] == valids_2[i]:
                        continue
                    valids1_str = f"{p1_object_list[object_idx_1]} and rank {valids_1[i]:d}"
                    valids2_str = f"{p1_object_list[object_idx_2]} and rank {valids_2[i]:d}"
                    if not valids2_str in remove_atoms:
                        if not valids1_str in remove_atoms:
                            remove_atoms.append(valids1_str)
        for ra in remove_atoms:
            pymol.cmd.remove(ra)
        create_list = list()
        for object_name in p1_object_list:
            if pymol.cmd.count_atoms(object_name) == 0:
                pymol.cmd.delete(object_name)
            else:
                create_list.append(object_name)
        pymol.cmd.create("input_p1", " or ".join(create_list))
        for object_name in create_list:
            pymol.cmd.delete(object_name)
        a, b, c, alpha, beta, gamma, spacegroup = pymol.cmd.get_symmetry("input")
        pymol.cmd.set_symmetry(
            "input_p1", 
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
            "P1")
        pymol.cmd.delete(f"input")
        for a in range(args.lattice_factors_x[0],args.lattice_factors_x[1]+1):
            for b in range(args.lattice_factors_y[0],args.lattice_factors_y[1]+1):
                for c in range(args.lattice_factors_z[0],args.lattice_factors_z[1]+1):
                    if verbose:
                        print("Generating unit cell",a,b,c)
                    object_list = symexpcell(
                        "i_",
                        "input_p1",
                        "input_p1",
                        args.duplicate_cutoff,
                        a, 
                        b, 
                        c)
                    for i in range(len(object_list)):
                        chain = string.ascii_uppercase[i]
                        chain = string.ascii_uppercase[i]
                        pymol.cmd.copy("mol%d" %object_count, object_list[i])
                        pymol.cmd.alter("mol%d" %object_count, 'chain="%s"' %chain)
                        pymol.cmd.group("global", "mol%d" %object_count)
                        if verbose:
                            pymol.cmd.save("mol%d_%d%d%d_%d_" %(input_idx,a,b,c,i) + args.output, "mol%d" %object_count)
                        pymol.cmd.delete(object_list[i])
                        object_count += 1
                    uc_count += 1
        pymol.cmd.delete("input_p1")

    pymol.cmd.save(args.output, "global")
