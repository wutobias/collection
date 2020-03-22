#!/bin/bash -f

######## Important executables ########
######## ~~~~~~~~~~~~~~~~~~~~~ ########

### Path to resp_it_mol2 scripts
RESP_IT_MOL2=/home/wulsdorf/shell-scripte/resp_it_mol2

resp_exe=$AMBERHOME/bin/resp

fconv_exe=/nfs/pub/file_utils/fconv

netcharge_exe=$RESP_IT_MOL2/calc_total_charge.py

### Can use different conformer generation strategy. Important thing is,
### that output file is multimol2 format with filename {input_filename}_conf.mol2
conf_gen_exe=$RESP_IT_MOL2/moe_stochastic_confs.sh

ensemble_filter_exe=$RESP_IT_MOL2/ensemble_filter_hb.py

GAUSS=/nfs/progs64/gaussian/g09c/g09/g09

ANTE=$AMBERHOME/bin/antechamber

ESPGEN=$AMBERHOME/bin/espgen

RESPGEN=$AMBERHOME/bin/respgen

CPPTRAJ=$AMBERHOME/bin/cpptraj

### This is specific for Klebe group computer systems.
### You might want to comment that.
source $RESP_IT_MOL2/login_g09.sh

########## F U N C T I O N S   F U N C T I O N S ##########
########## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##########

### very helpful function for setting and
### resetting gaussian scratch directory.
function scratch_set(){

    if [ ! -d scratch ]; then
        
        mkdir scratch
    
    fi
    
    export GAUSS_SCRDIR=$PWD/scratch

    
    #If $GAUSS_SCRDIR is present, empty it
    
    if [ ! $(find $GAUSS_SCRDIR -name "*" -type f | wc -l) -eq 0 ]; then
    
            rm -f $GAUSS_SCRDIR/*
    
    fi

}

### Working horse. Optimize geometry,
### calculate MEP and write to *gesp file

function calc(){

    qminput=$1
    
    qminput_name=${qminput%.*}

    netcharge=$2
    
    $ANTE -i $qminput -fi mol2 -o ${qminput_name}_opt.com -fo gcrt -ch ${qminput_name}_opt.chk -gn \%nproc\=$procs -gk \#HF/6-31G\*\ opt -nc $netcharge -pf y

    scratch_set

    $GAUSS ${qminput_name}_opt.com

    #extract final structure
    $ANTE -i ${qminput_name}_opt.log -fi gout -o tmp.mol2 -fo mol2 -nc $netcharge -pf y

    #get correct atomnames
    $ANTE -i tmp.mol2 -fi mol2 -o ${qminput_name}_opt.mol2 -fo mol2 -a $qminput -fa mol2 -ao name -nc $netcharge -pf y
    
    rm -f tmp.mol2  

    #write new gaussian file
    $ANTE -i ${qminput_name}_opt.log -fi gout -o ${qminput_name}_opt_esp.com -fo gcrt -ch ${qminput_name}_opt_esp.chk -gn \%nproc\=$procs -gk \#HF/6-31G\*\ pop\=mk\ iop\(6/33\=2\) -nc $netcharge -j 5 -pf y
            
    scratch_set

    $GAUSS  ${qminput_name}_opt_esp.com 
    
    $ESPGEN -i ${qminput_name}_opt_esp.log -o ${qminput_name}_opt.gesp
    
    scratch_set
    
}

################################# M A I N  S E C T I O N #################################
################################# ~~~~~~~~~~~~~~~~~~~~~~ #################################

confs=3
cutoff=1.
procs=2

### Amber and Gaussian checking
if [ ! $AMBERHOME ];then

    echo "Must set AMBERHOME."

    exit

fi

if ! type "g09" > /dev/null 2>&1; then

    echo "Gaussian not found."

    echo "Maybe source ${RESP_IT_MOL2}/login_g09.sh"

    exit
fi

### Command line checking...
if [ $# -lt 1 ]; then
    echo "Usage: resp_it.sh <mymolecule.mol2> NUMBEROFCONFORMERS RMSD-CUTOFF NUMBERPROCESSORS."
    exit
fi


if [ $# -eq 2 ]; then
    confs=$2
fi


if [ $# -eq 3 ]; then
    confs=$2
    cutoff=$3
fi

if [ $# -eq 4 ]; then
    confs=$2
    cutoff=$3
    procs=$4
fi

echo "Generating $confs conformers."
echo "Using cut-off=$cutoff Ang."
echo "Using $procs proc."


### read input

input=$1

input_name=${input%.*}

input_format=${input##*.}

if [ ! -d tmp_data ]; then

    mkdir tmp_data

fi
        
cd tmp_data

cp ../$input .

#Correct wrong atom types
#...in fact this can sometimes lead to wrong atom types?!
#$fconv_exe -T $input

#### Generate conformeres ####
#### ~~~~~~~~~~~~~~~~~~~~ ####

$conf_gen_exe $input

#First check, if we have found any additional conformers at all...
initial_confs=$(grep "@<TRIPOS>MOLECULE" ${input_name}_conf.mol2 | wc -l)

if [ $initial_confs -gt 1 ]; then

    #Fitler out conformeres with hbonds
    tmp=$RANDOM

    printf "parm ${input_name}_conf.mol2

        trajin ${input_name}_conf.mol2
        hbond out hbonds.dat
        go" > /tmp/${tmp}_hbonds.in
    
    $CPPTRAJ -i /tmp/${tmp}_hbonds.in

    if [ -f hbonds.dat ]; then
        grep -v "#" hbonds.dat > tmp.dat
        mv -f tmp.dat hbonds.dat
        awk '{print $1-1 " " $2}' hbonds.dat > tmp.dat
        mv -f tmp.dat hbonds.dat
    else
        touch hbonds.dat
    fi
   
    $ensemble_filter_exe hbonds.dat ${input_name}_conf.mol2
    
    #Cluster residual conformers
    $fconv_exe -clust ${input_name}_conf_nohb.mol2 --t=${input_name}_conf_clustered.mol2 --r=$cutoff
    
        #It is possible that less than $confs have been found, so determine number of actually found conformers
        total_confs=$(grep "@<TRIPOS>MOLECULE" ${input_name}_conf_clustered.mol2 | wc -l)
        
    #If we found less conformations than we intended to, use all of them.
        if [ $total_confs -lt $confs ]; then
        
            confs=$total_confs
            
    fi

    #Extract the conformations from the clustering
    for (( k=0; k<$confs; k++ )); do

        $fconv_exe --e=$k ${input_name}_conf_clustered.mol2 --t=${input_name}_conf${k}.mol2
    
    done
        
else

    confs=1
    ln -s ${input_name}_conf.mol2 ${input_name}_conf0.mol2

fi


if [ -f ${input_name}.gesp ]; then

    rm -f ${input_name}.gesp
    
fi


##### Calculate Conformers via QM #####
##### ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #####

netcharge=$($netcharge_exe ${input_name}.mol2)

for (( k=0; k<$confs; k++ ))
do

    calc ${input_name}_conf${k}.mol2 $netcharge

    cat ${input_name}_conf${k}_opt.gesp >> ${input_name}.gesp

done

##### Resp fit all MEPs #####
##### ~~~~~~~~~~~~~~~~~ #####

#Go!
$ANTE -i ${input_name}.mol2 -fi mol2 -o ${input_name}.ac -fo ac -nc $netcharge -pf y 

$RESPGEN -i ${input_name}.ac -o ${input_name}.resp1 -f resp1 -n $confs

$RESPGEN -i ${input_name}.ac -o ${input_name}.resp2 -f resp2 -n $confs

$resp_exe -O -i ${input_name}.resp1 -o ${input_name}.resp1out -e ${input_name}.gesp -t ${input_name}_resp1.crg 

$resp_exe -O -i ${input_name}.resp2 -o ${input_name}.resp2out -e ${input_name}.gesp -t ${input_name}_resp2.crg -q ${input_name}_resp1.crg 

$ANTE -i ${input} -fi mol2 -o ${input_name}_resp.mol2 -fo mol2 -cf ${input_name}_resp2.crg -c rc -nc $netcharge -at sybyl -pf y

cp ${input_name}_resp.mol2 ..

for (( k=0; k<$confs; k++ ))
do

    $ANTE -i ${input_name}_conf${k}_opt.mol2 -fi mol2 -o ${input_name}_conf${k}_resp.mol2 -fo mol2 -cf ${input_name}_resp2.crg -c rc -nc $netcharge -at sybyl -pf y

    cp ${input_name}_conf${k}_resp.mol2 ..

done

cd ..
tar cjf ${input_name}_tmp_data.tar.bz tmp_data
### This is dangerous, uncomment only when neccesary.
#rm -rf tmp_data