#!/bin/bash

if [ ! $# -eq 5 ]; then

	echo "Usage: compare_halfs.sh <PARM> <REF> <TRAJ1> <TRAJ2> <MASK>"
	exit 0

fi

if [ ! $CPPTRAJ ]; then

        echo "Setting CPPTRAJ\=\"cpptraj.OMP\" "
        CPPTRAJ="cpptraj.OMP"
        
else

        echo "CPPTRAJ is set to $CPPTRAJ"

fi

prmtop=$1
ref=$2
traj1=$3
traj2=$4
mask=$5

printf "
parm $prmtop
reference $ref $mask [ref]
trajin $traj1
strip !($mask)
rmsd ref [ref] rmsd_1 
go

clear trajin

trajin $traj2
strip !($mask)
rmsd ref [ref] rmsd_2
go

clear trajin

trajin $traj1
trajin $traj2
strip !($mask)
rmsd ref [ref] rmsd_total out rmsd_series.dat
go

clear trajin

kde rmsd_1     out rmsd_1_kde.dat     name kde_1 step 0.05
kde rmsd_2     out rmsd_2_kde.dat     name kde_2 step 0.05
kde rmsd_total out rmsd_total_kde.dat name kde_t step 0.05
go
divergence ds1 kde_1 ds2 kde_2
go
" > compare_halfs.cpptraj

$CPPTRAJ -i compare_halfs.cpptraj
/home/wulsdorf/python-scripte//get_kb.py rmsd_series.dat