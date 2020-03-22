#!/usr/bin/env python

import sys
import os

fconv = "/nfs/pub/file_utils/fconv"

if len(sys.argv) != 3:
    print "Usage: ensemble_filter_hb.py <ID_vs_hydrogenbondcounts.dat> <ensemble.mol2>"
    exit()

for f in sys.argv:
    if not os.path.exists(f):
        print "%s not found. Exit." %f

### Hbond information for ensemble from hbond_information file as obtained from vmd or cpptraj (frame vs. number of hbonds)
### 0 0
### 1 0
### 2 0
### 3 0
### 4 0
### 5 0
### 6 0
### 7 2
### 8 0
### 9 2
### 10 0
### 11 0

hbond_counts_tmp  = open(sys.argv[1],"r")
hbond_counts_file = hbond_counts_tmp.readlines()
hbond_counts_tmp.close()

outname = sys.argv[2].replace(".mol2", '') + "_nohb.mol2"

ensemble_count = 0
hbond_limit    = 0

while ensemble_count < 1:

    n_lines = 0
    for line in hbond_counts_file:

        if len(line.rstrip().split()) == 0:
            continue
        
        if line[0] == '#':
            continue
        
        n_lines += 1

        extract_in = ''
        merge_in   = ''

        if int(line.rstrip().split()[1]) == hbond_limit:

            #extract it
            extract_in += " --e=" + str(int(line.rstrip().split()[0]))
            extract_in += " " + sys.argv[2]
            extract_in += " --t=tmp_mol.mol2"
            os.system (fconv + extract_in)

            if ensemble_count == 0:
        
                os.system("mv tmp_mol.mol2 " + outname)
                ensemble_count += 1
                continue

            #merge it
            merge_in   += " -m "  + outname + " tmp_mol.mol2"
            merge_in   += " --t=tmp_ensemble.mol2"
            os.system (fconv + merge_in)
            os.system ("cp -f tmp_ensemble.mol2 " + outname)
        
            ensemble_count += 1

    if n_lines==0:
        break
        
    hbond_limit += 1
    
#clean
if n_lines>0:
    os.system ("rm -f tmp_ensemble.mol2 tmp_mol.mol2")
else:
    os.system("cp -f %s %s" %(sys.argv[2], outname))