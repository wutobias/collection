#!/usr/bin/env python

import sys

path2molecule = sys.argv[1]

file_reference = open(path2molecule,"r")

reference = file_reference.readlines()

file_reference.close()

#'index version' to get the block between @<TRIPOS>ATOM and @<TRIPOS>BOND

for index,line in enumerate(reference):

	if line == "@<TRIPOS>ATOM\n":

		start_ATOM = index+1

	if line == "@<TRIPOS>BOND\n":

		end_ATOM = index

		break

total_charge = 0.0 

for i, line in enumerate(reference[start_ATOM:end_ATOM]):

	c = line.rstrip().split()

	if c[-1].endswith("-"):

		total_charge -= float(c[-1].replace('-',''))

	else:

		total_charge += float(c[-1].replace('+',''))

print int(round(total_charge))