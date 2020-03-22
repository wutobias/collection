#!/usr/bin/env python

import sys

if len(sys.argv) < 2:

        print "Usage: esptotbl.py my_esp_fit.gesp > my_esp_fit.dat"
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

print 'Xcrd \t Ycrd \t Zcrd \t pot'
for key in esp_dict.keys():
  print "%6.3f \t %6.3f \t %6.3f \t %6.3f" %(key[0], key[1], key[2], esp_dict[key])
print "\n"