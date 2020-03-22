#!/usr/bin/env python

import sys

if len(sys.argv) < 2:

        print "Usage: esptopdb.py my_esp_fit.gesp > my_esp_fit.pdb"
        exit(0)

esp_file = open(sys.argv[1], "r")
esp      = esp_file.readlines()
esp_file.close()

# Bohr radius in 10^-10 m
scale = 0.5291772

input_name=sys.argv[1]

number_atoms = int(esp[0].rstrip().split()[0])
esp_dict     = {}
mol_dict     = {}

for i, line in enumerate( ( esp[ 1 : ( number_atoms ) ] ) ):

        c = line.rstrip().split()
        mol_dict[ float(c[0]) * scale, float(c[1]) * scale, float(c[2]) * scale ] = 0.00

for i, line in enumerate( ( esp[ ( number_atoms+1 ) : ] ) ):

        c = line.rstrip().split()
        esp_dict[  float(c[1]) * scale, float(c[2]) * scale, float(c[3]) * scale ] = float(c[0]) * 10**4

print 'REMARK default_name'

#for key in mol_dict.keys():
#  print '%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s' %('HETATM',i+1,'O','', 'MOL', 'A', i+1, '', key[0]*scale, key[1]*scale, key[2]*scale, 0.00, mol_dict[key], 'O', '')

for key in esp_dict.keys():
  print '%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s' %('HETATM',i+1,'O','', 'HOH', 'B', i+1, '', key[0], key[1], key[2], 0.00, esp_dict[key], 'O', '')

print "END"