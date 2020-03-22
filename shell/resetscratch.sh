#!/bin/bash

R=$RANDOM

if [ ! -d $PWD/scratch ]; then 
  mkdir /scratch/scratch${R}
fi

export GAUSS_SCRDIR=/scratch/scratch${R}
