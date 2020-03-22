#!/bin/bash -f

input=$1
input_name=${input%.*}
smi=$(cat $input)

moebatch_exe=/nfs/progs64/ccg/moe64/bin/moebatch

$moebatch_exe -exec "\
mol = sm_BuildMOL '$smi'; \
mol_Create mol ; \
issues = StructurePreparation [cmd: 'cli']; \
atoms = Atoms[]; \
pot_Load  '/nfs/progs/ccg/moe/lib/mmff94x.ff'; \
Protonate3D [atoms,atoms,atoms,[],[],[], verbose: 0]; \
MM [ pot_charge:1 ]; \
WriteTriposMOL2 '$input_name.mol2'; \
Close []; \
"
