#!/bin/bash -f

if [ $# -lt 1 ]; then

	echo "Usage: count_wat.sh protein1.pdb protein2.pdb ..."
	exit 1

fi

PDB_LIST=$*
PDB_LIST=${PDB_LIST#count_wat.sh}

for PDB in $PDB_LIST; do

	C=$(grep -H 'WAT\|HOH' $PDB | grep 'ATOM\|HETATM' | grep "O  " | wc -l)
	echo "$PDB: $C waters"

done