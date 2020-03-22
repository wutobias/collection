#!/usr/bin/env python

import argparse
import numpy as np
import string
import __main__
import xtalmd
import os
from xtalmd.utils import cellbasis

parser = argparse.ArgumentParser(description="Rotate structures according to a reference rotation.\
                                              Written by Tobias Wulsdorf, AG Klebe, tobias.wulsdorf@gmail.com")

parser.add_argument('--reference',       '-r',  type=str,\
                                                help='Input PDB file used as target structure for superposition.',\
                                                required=True)
parser.add_argument('--input',           '-i',  type=str,\
                                                help='Reference PDB file. Is superimposed with input pdb structure \
and then used as template of object matrix and transformation.',\
                                                required=True)
parser.add_argument('--mobile',          '-m',  type=str,\
                                                help='Mobile structures which will be aligned according to object \
matrix file in --reference.',\
                                                required=True,
                                                nargs='+')
parser.add_argument('--prefix',          '-pre',  type=str,\
                                                help='Prefix for output files',\
                                                default='align')

if __name__ == "__main__":

    args = parser.parse_args()

    __main__.pymol_argv = ['pymol','-c'] # Pymol: quiet and no GUI
    #__main__.pymol_argv = ['pymol','-qQc'] # Pymol: quiet and no GUI
    import pymol
    pymol.finish_launching()

    pymol.cmd.set("pdb_conect_all", "on")
    pymol.cmd.set("connect_mode", "1")
    pymol.cmd.set("retain_order", "on")
    pymol.cmd.set("pdb_retain_ids", "on")

    pymol.cmd.load(args.reference, "reference")
    pymol.cmd.load(args.input, "input")

    pymol.cmd.do("align input and name CA, reference and name CA")

    for i, mobile in enumerate(args.mobile):
    
        base = os.path.basename(mobile)

        pymol.cmd.load(mobile, "strc_%d" %i)
        pymol.cmd.matrix_copy("input", "strc_%d" %i)
        pymol.cmd.save("%s."%args.prefix + base, "strc_%d" %i)