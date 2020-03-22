#!/bin/bash -f

if [ $# -lt 1 ]; then
	echo "Usage: obabel_conformers.sh <mymolecule.mol2>"
	exit
fi



if [ $# -lt 2 ]; then
	echo "Generating 3 conformers."
	confs=3
fi



if [ $# -eq 2 ]; then
	echo "Generating $2 conformers."
	confs=$2
fi



i=$1
j=${i%.*}

obabel $i -O ${j}_conf.mol2 --conformer --nconf $confs --score energy --writeconformers -h
