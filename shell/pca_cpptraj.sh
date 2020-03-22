#!/bin/bash -f

###
### T. Wulsdorf
### AG Klebe, Marburg, 2016
### tobias.wulsdorf@uni-marburg.de
###

if [ ! $# -eq 4 ]; then

	echo "pca_cpptraj.sh [prmtop] [traj] [mask] [number-of-evecs]"
	exit 1

fi

cpptraj="cpptraj"
PRMTOP=$1
TRAJ=$2
MASK=$3
EVECS=$4

printf """
#Step 1: Calculation of the coordinate covariance matrix
#-------------------------------------------------------

strip !$MASK
autoimage
rms first
trajout prod.nc
go
clear all
parm $PRMTOP
parmstrip !$MASK
parmwrite out prmtop.prmtop
go
clear all
parm prmtop.prmtop
trajin prod.nc
average crdset mol-avg
createcrd mol-traj
run
crdaction mol-traj rms ref mol-avg
crdaction mol-traj matrix covar name mol-covar

#Step 2: Calculate principal components and coordinate projections
#-----------------------------------------------------------------

runanalysis diagmatrix mol-covar out mol-evecs.dat \
    vecs $EVECS name myEvecs \
    nmwiz nmwizvecs $EVECS nmwizfile mol.nmd
crdaction mol-traj projection mol-proj modes myEvecs \
  beg 1 end $EVECS crdframes 0,last out proj.dat

#Step 3: Visualizing principal components
#----------------------------------------

clear all
readdata mol-evecs.dat name Evecs
parm $PRMTOP
parmstrip !($MASK)
parmwrite out mol-modes.prmtop

runanalysis modes name Evecs eigenval beg 1 end $EVECS out evecs-loadings.dat
runanalysis modes name Evecs fluct    beg 1 end $EVECS out evecs-fluct.dat

""" > pca.cpptraj

for ((i=0;i<$EVECS;i++)); do

  a=$(python -c "print ${i}+1")
  printf """
runanalysis modes name Evecs trajout mol-mode${a}.nc \
pcmin -100 pcmax 100 tmode ${a} trajoutfmt netcdf\
""" >> pca.cpptraj

done

$cpptraj -p $PRMTOP -y $TRAJ -i pca.cpptraj

exit 0
