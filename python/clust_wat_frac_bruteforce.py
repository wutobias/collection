#!/usr/bin/env python

### Written by T. Wulsdorf @ AG Klebe Marburg University
### Please contact via tobias.wulsdorf@uni-marburg.de

### Use field classes abd grid_water
from io_operations import load_dx, write_files
from box_crds import box_crds 

### Other
import argparse
import multiprocessing
import matplotlib as mpl
import time
import os
import mdtraj as md

### Plotting
mpl.use('Agg')
import matplotlib.pyplot as plt

### Numerical
import numpy as np
import scipy.spatial as spatial_scipy

parser = argparse.ArgumentParser(description='Apply brute force WaterMap clustering to find hydration sites. Written by Tobias Wulsdorf, \
	                                          tobias.wulsdorf@gmail.com')

parser.add_argument('--traj',            '-t',   type=str,   help='Trajectory file', required=True)
parser.add_argument('--topo',            '-p',   type=str,   help='Parameter file in pdb format', required=True)
parser.add_argument('--strc_ref',        '-tr',  type=str,   help='Refrence structure file used for rms fitting prior to HSA.', required=False)
parser.add_argument('--topo_ref',        '-pr',  type=str,   help='Parameter file for rms fitting.', required=False)
parser.add_argument('--prefix',          '-pf',  type=str,   help='Filename prefix', default="HSA")
parser.add_argument('--start',           '-s0',  type=int,   help='First frame', required=True)
parser.add_argument('--stop',            '-s1',  type=int,   help='Last frame',  required=True)
parser.add_argument('--stride',          '-st',  type=int,   help='Stride frames', default=100)
parser.add_argument('--center',          '-c',   type=float, nargs=3, help='Grid box center.', default=[0.,0.,0.])
parser.add_argument('--dim',             '-d',   type=float, nargs=3, help='Grid box dimensions in angstrom', default=[20.,20.,20.])
parser.add_argument('--radius',          '-r',   type=float, help='Cluster center radius.', default=1.0)
parser.add_argument('--threshold',       '-th',  type=float, help='Threshold for algorithm to stop. The value of threshold must be given in multiples of bulk density. Default 2.0', default=2.0)
parser.add_argument('--nproc',           '-np',  type=int,   help='Number of cores used for calculation.', default=2)
parser.add_argument('--grid_files',      '-gf',  type=str,   help='Directory containing center and frac2real information for each frame.')
parser.add_argument('--reference_frame', '-rf',  type=int,   help='Reference frame used for real coordinate output in case of rotation grids.', default=0)


def main():

  out_q = multiprocessing.Queue()
  processes = list()
  # Dummy crds array
  crds      = np.zeros((1,3))
  cluster_centers = list()

  external_grid = False
  center        = np.zeros(3)
  dim           = np.zeros(3)

  fitting = False

  if args.strc_ref != None or args.topo_ref != None:
    fitting = True
    if args.strc_ref == None or args.topo_ref == None:
      print "--strc-ref and --topo-ref must be given both."
      print "Exit."
      quit()

  #Adjusted Bulk density is total number of water molecules
  #expected in one hydration site in bulk for analyzed number 
  #of frames in this run.
  adjust_bulk_density = 0.0334 * (4./3. * np.pi * args.radius**3) * float((args.stop - args.start)/args.stride)
  print "Adjusted number of waters per HS in bulk: ", adjust_bulk_density

  if args.grid_files != None:
    external_grid = True

  elif args.center == None:
    raise Error("If external grid alignment not used, must specify grid center.")

  else:
    dim = np.array(args.dim)

  def worker(out_q, _range):

    g = box_crds(args.topo, np.array(args.dim))
    c = list()

    with md.open(args.traj) as md_traj:
      for i in _range:

        md_traj.seek(i)
        frame   = md_traj.read_as_traj(g.topo, n_frames=1, stride=1)

        if external_grid:

          ### Unit is cell matrix vectors. It will be transformed
          ### into frac2real within the set_frame call. Therefore we
          ### MUST assume that cell matrix vectors are perfectly normalized.
          if os.path.exists(args.grid_files+"/%d_unit.dat.gz"%i):
            unit  = np.loadtxt(args.grid_files+"/%d_unit.dat.gz"%i)
          else:
            unit  = np.loadtxt(args.grid_files+"/%d_unit.dat"%i)

          if os.path.exists(args.grid_files+"/%d_center.dat.gz"%i):
            center = np.loadtxt(args.grid_files+"/%d_center.dat.gz"%i)
          else:
            center = np.loadtxt(args.grid_files+"/%d_center.dat"%i)

          g.set_frame(frame, unit, center)

        else:
          g.set_frame(frame, np.eye(3,3), args.center)

        c.append(g.get_inside_crds_frac())

      out_q.put(c)

  _range           = range(args.start, args.stop, args.stride)
  _start           = 0
  _stop            = 0
  _steps_p_proc    = int(len(_range)/args.nproc)
  _addsteps_p_proc = int(len(_range)%args.nproc)

  print "Processing with %s cores..." %args.nproc

  for p in range(args.nproc):

    _start = _stop
    _stop  = _start + _steps_p_proc

    if _addsteps_p_proc > 0:

      _stop            += 1
      _addsteps_p_proc -= 1

    p = multiprocessing.Process(target=worker, args=(out_q, _range[_start:_stop]))
    processes.append(p)
    p.start()

  non_empty_found = False
  for p in processes:

    c = out_q.get()
    if len(c) != 0:
      if not non_empty_found:
        crds = c[0]
        for i in range(1,len(c)):
          crds_new = np.append(crds, c[i], axis=0)
          crds     = crds_new
        non_empty_found = True
      else:
        for i in range(0,len(c)):
          crds_new = np.append(crds, c[i], axis=0)
          crds     = crds_new

  for p in processes:

    p.join()

  print "Clustering..."
  L                         = crds.shape[0]                           # Holds total number of water coordinates, max(i)
  cluster_centers_idx       = np.array([], dtype=int)                 # Holds index of cluster center coordinates
  cluster_centers_idx_count = np.array([], dtype=float)               # Holds total number of coordinates in each cluster center

  Tree               = spatial_scipy.cKDTree(crds)                    # KDTree used for distance calculation
  NN_list            = Tree.query_ball_tree(Tree, 2.*args.radius)        # Nearest Neighbour list constructred from KDTree. It is a list A of lists, 
                                                                      # where each coordinate index A[i] has a list that holds the indices B[j]
                                                                      # of its nearest neighbours within args.radius.
  NN_count           = np.zeros(L, dtype=int)                         # Total number of neighbours found within args.radius of each coordinate i
  valids             = np.ones(L, dtype=bool)                         # List of length max(i) holding bool(False) at each coordinate that has not been excluded

  for i, NN in enumerate(NN_list):

    NN_count[i] = len(NN)

  max_i   = np.argmax(NN_count)
  max_i_N = float(len(NN_list[max_i]))

  while max_i_N > args.threshold * adjust_bulk_density:

    cluster_centers_idx_tmp = np.append(cluster_centers_idx, max_i)
    cluster_centers_idx     = cluster_centers_idx_tmp

    cluster_centers_idx_count_tmp = np.append(cluster_centers_idx_count, max_i_N)
    cluster_centers_idx_count     = cluster_centers_idx_count_tmp

    #print "New cluster center is ", max_i
    #print "...with ", max_i_N, " waters."

    valids[NN_list[max_i]] = False
    valids[max_i]          = False
    valids_idx             = np.where(valids)[0]
    N_valids               = valids_idx.shape[0]
    if N_valids == 0:
      #print "Did not find any new cluster center."
      #print "Exit."
      break

    max_i            = np.argmax(NN_count[valids])
    max_i            = valids_idx[max_i]
    max_i_N          = float(NN_count[max_i])

  cluster_centers = crds[cluster_centers_idx]

  print "Writing files..."
  remarks="Clustered with clust_wat_frac_watermap.py. "
  remarks+="Local current time: %s. " %time.asctime( time.localtime(time.time()) )
  remarks+="Clustersize=%s " %args.radius
  remarks+="traj=%s " %args.traj
  remarks+="topo=%s " %args.topo
  remarks+="traj=%s " %args.traj
  remarks+="prefix=%s " %args.prefix
  remarks+="start=%s " %args.start
  remarks+="stop=%s " %args.stop
  remarks+="stride=%s " %args.stride
  if args.grid_files != None:
    remarks+="grid_files=%s " %args.grid_files
  else:
    remarks+="center=%s " %args.center
    remarks+="dim=%s " %args.dim

  args.reference_frame = 0
  g = box_crds(args.topo, np.array(args.dim), args.strc_ref, args.topo_ref)
  frame = md.load_frame(args.traj, args.reference_frame, top=g.topo)
  
  if external_grid:

    if os.path.exists(args.grid_files+"/%d_unit.dat.gz"%i):
      unit  = np.loadtxt(args.grid_files+"/%d_unit.dat.gz"%i)
    else:
      unit  = np.loadtxt(args.grid_files+"/%d_unit.dat"%i)

    if os.path.exists(args.grid_files+"/%d_center.dat.gz"%i):
      center = np.loadtxt(args.grid_files+"/%d_center.dat.gz"%i)
    else:
      center = np.loadtxt(args.grid_files+"/%d_center.dat"%i)

    g.set_frame(frame, unit, center)

  else:
    g.set_frame(frame, np.eye(3,3), args.center)
      
  write_files(XYZ=cluster_centers/g.f.bins,               
              Format='PDB', 
              Filename=args.prefix+"_cluster-centers-frac.pdb", 
              Remarks=remarks, 
              Value=cluster_centers_idx_count/adjust_bulk_density)
  write_files(XYZ=g.f.get_real(cluster_centers), 
              Format='PDB', 
              Filename=args.prefix+"_cluster-centers.pdb",      
              Remarks=remarks, 
              Value=cluster_centers_idx_count/adjust_bulk_density)
  #Debug
  #write_files(XYZ=g.f.get_real(crds),            Format='PDB', Filename="raw_data.pdb",                     Remarks=remarks)
  g.frame.save_pdb(filename=args.prefix+'_reference_frame_%d.pdb' %args.reference_frame, force_overwrite=True)

  print "Done."
  
if __name__ == "__main__":

  args = parser.parse_args()

  main()