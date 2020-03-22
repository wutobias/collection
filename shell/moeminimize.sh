#!/bin/bash -f

ligandname=${1%.*}

/nfs/progs64/ccg/moe64/bin/moebatch -exec "\
ReadTriposMOL2 '${1}'; \
pot_Load  '/nfs/progs/ccg/moe/lib/mmff94x.ff';
MM [ pot_charge:0, keep_chirality:'geometry' ];
WriteTriposMOL2 '${ligandname}_min.mol2'; \
Close []; "
