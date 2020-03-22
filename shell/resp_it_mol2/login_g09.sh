#!/bin/bash -f

if [ $USER!=root ]; then
   export wd=$PWD
   export g09root=/nfs/progs64/gaussian/g09c
   if [ ! -d ${HOME}/.gaussian ]; then
      mkdir $HOME/.gaussian
      mkdir $HOME/.gaussian/scratch
   fi

   k=$(uname -r)
   if [[ $k == 3.2.* ]]; then
      if [ "$LD_LIBRARY_PATH" ]; then
         export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/nfs/progs64/libc-2.15
      else
         export LD_LIBRARY_PATH=/nfs/progs64/libc-2.15
      fi
   fi

   export GAUSS_SCRDIR=$HOME/.gaussian/scratch   
   cd $g09root/g09/bsd
   source $g09root/g09/bsd/g09.profile
   cd $wd
fi
