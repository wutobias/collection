#!/bin/bash -f

i=$1
j=${i%.*}

moebatch=/nfs/progs64/ccg/moe64/bin/moebatch

$moebatch -exec "\
    ReadPDB '${j}.pdb'; 
    issues = StructurePreparation [cmd: 'cli', reportName:'${j}_moe.log', batch_protonate3d: 1]; 
    atoms = Atoms[]; 
    Protonate3D [atoms,atoms,atoms,[],[],[], verbose: 0]; 
    WritePDB ['${j}_moe.pdb']; 
    Close []; 
    "
fconv -W ${j}_moe.pdb --t=${j}_nowat.pdb
fconv -w ${j}_moe.pdb --t=${j}_wat.mol2

fconv -l ${j}_nowat.pdb
fconv -L ${j}_nowat.pdb --t=${j}_nowat_prot.pdb
renumberpdb.py ${j}_nowat_prot.pdb > ${j}_nowat_prot_re.pdb
amberready.sh ${j}_nowat_prot_re.pdb
if [ -f ${j}_nowat_prot_re_amb.pdb ]; then
	ln -s ${j}_nowat_prot_re_amb.pdb prot.pdb
fi
