#!/bin/bash -f

i=$1
j=${i%.*}

moebatch=/nfs/progs64/ccg/moe64/bin/moebatch

$moebatch -exec "\
    ReadTriposMOL2 '$i'; \
    ConfSearch [method: 'Stochastic', outfile: '${j}.mdb']; \
    db_ExportTriposMOL2 ['${j}.mdb', '${j}_conf.mol2']; \
    Close []; "
    
#WriteTriposMOL2 ['${j}_moe.mol2', ]; \ 
