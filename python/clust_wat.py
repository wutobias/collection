#!/usr/bin/env python

### Written by T. Wulsdorf @ AG Klebe Marburg University
### Please contact via tobias.wulsdorf@uni-marburg.de

### Numerical libraries
import numpy as np
from scipy import spatial as spatial_scipy

### Use field classes abd grid_water
from io_operations import load_dx, write_files

### Other
import sys
import time

if len(sys.argv) != 5:

  print "Usage: clust_wat.py <field.dx> <threshold> <size> <cluster-centers.pdb>"
  exit()

file        = sys.argv[1]
threshold   = float(sys.argv[2])
clustersize = float(sys.argv[3])
output      = sys.argv[4]

print "Loading file ..."
f,v    = load_dx(file)
valids = np.where(v >= threshold)
cluster_centers = list()

print "Clustering ..."
while valids[0].shape[0] > 0:

  local_max_x = np.where(v == np.amax(v))
  local_max_x = f.get_real(np.array([ local_max_x[0][0],
                                      local_max_x[1][0],
                                      local_max_x[2][0]]))

  ### Reshape valids array such that
  ### [[0,0,0],
  ###  [1,1,0]
  ###  [0,1,0],
  ###  ...,
  ###  ]
  valids_reshape = np.stack([valids[0], valids[1], valids[2]], axis=1)
  invalids       = np.where(spatial_scipy.distance.cdist(f.get_real(valids_reshape), np.array([local_max_x])) <= clustersize)[0]
  v[valids_reshape[invalids][:,0], valids_reshape[invalids][:,1], valids_reshape[invalids][:,2]] = 0.

  cluster_centers.append(local_max_x)

  valids = np.where(v >= threshold)

if len(cluster_centers) == 0:

  print "No cluster centers found."
  exit(0)

print "Writing cluster center pdb..."
remarks="Clustered with clust_wat.py. "
remarks+="Local current time: %s. " %time.asctime( time.localtime(time.time()) )
remarks+="Threshold=%s " %threshold
remarks+="Clustersize=%s"   %clustersize
write_files(XYZ=np.asarray(cluster_centers), Format='PDB', Filename=output, Remarks=remarks)
print "Done."
