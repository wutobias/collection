#!/usr/bin/env python

from watcor import watcor2
import sys

if len(sys.argv) != 4:
   print "Usage: watcor_analysis.py <TRAJ> <PARM> <hydration-sites>"
   exit(1)

traj=sys.argv[1]
parm=sys.argv[2]
hs=sys.argv[3]

wc2=watcor2.watcor(traj, parm, hs)
corrmat=wc2.cor_pop(radius=1.5)
wc2.pymol_out(corrmat)
