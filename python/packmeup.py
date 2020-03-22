#!/usr/bin/env python

### Written by Tobias Wulsdorf, AG Klebe, 2017
### Please contact via tobias.wulsdorf@staff.uni-marburg.de

import pytraj as pt
import numpy as np
from scipy.spatial import distance
from pygist import io
import argparse

from mystats import bootstrap_resample

parser = argparse.ArgumentParser(description='Extract ligands from ensemble/trajectory and place it randomly around a protein.')

parser.add_argument('--traj',        '-t',   type=str,   help='Ligand Trajectory file'                          , required=True)
parser.add_argument('--prmtop',      '-p',   type=str,   help='Ligand Parameter file'                           , required=True)
parser.add_argument('--prefix',      '-pf',  type=str,   help='Prefix name for output'                          , default="Molecule")
parser.add_argument('--outformat',   '-of',  type=str,   help='Output format Molecules.'                        , default='pdb')
parser.add_argument('--prot_strc',   '-ps',  type=str,   help='Protein structure.'                              , required=True)
parser.add_argument('--prot_prmtop', '-pp',  type=str,   help='Protein Parameter.'                              , required=True)
parser.add_argument('--num_mol',     '-n',   type=int,   help='Number of molecules to generate'                 , default=10)
parser.add_argument('--buff_pro',    '-bp',  type=int,   help='Buffer distance protein'                         , default=4.0)
parser.add_argument('--buff_lig',    '-bl',  type=int,   help='Buffer distance ligand'                          , default=3.0)
parser.add_argument('--debug',       '-db',  type=int,   help='Writes some additional files.'                   , default=0)
parser.add_argument('--sample_rate', '-sp',  type=int,   help='Sampling rate for phi/theta and rotation angles.', default=360)

DEG2RAD        = np.pi/180.
RAD2DEG        = 180./np.pi

def main():

	print "Precomputing some stuff..."

	### Some general parameters
	N_ligs         = args.num_mol
	_buffer_pro    = args.buff_pro
	_buffer_lig    = args.buff_lig
	angle_sampling = args.sample_rate
	iteration      = 0
	debug          = bool(args.debug)
	prefix         = args.prefix
	
	### Load files
	lig         = pt.load(args.traj, args.prmtop)
	prot        = pt.load(args.prot_strc, args.prot_prmtop)
	
	### Bootstrap ligand conformations
	frame_array = bootstrap_resample(np.arange(lig.n_frames), N_ligs)
	
	### Angle array. We use it for spherical coordinate transformations
	### and rotation matrices
	theta   = np.linspace(0.0, 180.0, angle_sampling, endpoint=False)
	phi     = np.linspace(-180.0, 180.0, angle_sampling, endpoint=False)
	
	### Precompute coordinates for ligands, protein and whole box
	xyz_ligand  = lig[":1"].xyz[frame_array]
	xyz_prot    = prot["!(:WAT,Na+,Cl-)"].xyz
	xyz_box     = prot.xyz
	
	### Precompute geometry centers for ligands, protein and whole box
	cog_pro = xyz_prot.mean(axis=1)
	cog_box = xyz_box.mean(axis=1)
	cog_lig = xyz_ligand.mean(axis=1)
	cog_lig = np.expand_dims(cog_lig, axis=1)
	
	### Precompute largest sphere that could still encompass protein,
	### box or ligand
	r_pro  = np.max(np.linalg.norm(cog_pro-xyz_prot, axis=2))
	r_box  = np.max(np.linalg.norm(cog_box-xyz_box, axis=2))
	r_lig  = np.max(np.linalg.norm(cog_lig-xyz_ligand, axis=2), axis=1)
	
	### sanity check
	#if r_box < np.max(r_pro + _buffer_pro + 2.*r_lig):
	#    raise ValueError("Box is too small. Make bigger box or decrease buffer.")
	
	### dist_lig: distance between every ligand and its nearest neighbour
	### min_dist: sum of the radii of every ligand and its nearest neighbour
	dist_lig =    r_lig
	min_dist = 2.*r_lig
	
	### Translate ligands to origin
	xyz_ligand_origin = xyz_ligand - cog_lig

	print "Finding optimal ligand positions..."
	    
	### Assume each ligand is a sphere with radius r_lig
	### and place it such that dist_lig > min_dist
	while not np.all(dist_lig > min_dist):
	    
	    ### Generate random positions on sphere
	    theta_random_idx = np.random.randint(angle_sampling, size=N_ligs)
	    phi_random_idx   = np.random.randint(angle_sampling, size=N_ligs)
	    
	    theta_random     = theta[theta_random_idx] * DEG2RAD
	    phi_random       = phi[phi_random_idx]     * DEG2RAD
	    
	    rad              = r_pro + _buffer_pro + r_lig
	    
	    x_sphere, y_sphere, z_sphere = rad*np.sin(theta_random)*np.cos(phi_random),\
	                                   rad*np.sin(theta_random)*np.sin(phi_random),\
	                                   rad*np.cos(theta_random)
	            
	    xyz_sphere         = np.stack((x_sphere,y_sphere,z_sphere), axis=1)
	    
	    ### Translate sphere centers to protein cog
	    xyz_sphere         = xyz_sphere + cog_pro
	    # Debug
	    if debug:
	    	io.write_files(XYZ=xyz_sphere, Format="PDB", Filename=prefix+"_"+str(iteration)+"_principal_positions_centers.PDB")
	    xyz_sphere         = np.expand_dims(xyz_sphere, axis=1)
	
	    xyz_ligand_new     = xyz_sphere + xyz_ligand_origin
	    
	    # Debug
	    if debug:
	    	pt.write_trajectory(prefix+"_"+str(iteration)+"_principal_positions_strct.pdb", xyz_ligand_new, top=lig.top, overwrite=True)
	
	    cog_ligand_new     = xyz_ligand_new.mean(axis=1)
	    
	    ### Final distance check
	    dists_lig         = distance.cdist(cog_ligand_new, cog_ligand_new)
	    np.fill_diagonal(dists_lig, np.inf)
	    dist_lig_arg      = np.argmin(dists_lig, axis=1)
	    dist_lig          = dists_lig[(np.arange(N_ligs), dist_lig_arg)]
	    min_dist          = r_lig + r_lig[dist_lig_arg] + _buffer_lig*0.5
	    
	    if debug:

	    	print "Iteration: ", iteration
	    	print "Distances: ", dist_lig
	    	print "Min Dists: ", min_dist

	    iteration += 1
	    
	print "Generating random orientations..."

	### Give random orientation
	Rx               = np.zeros((N_ligs, 3, 3), dtype=float)
	Ry               = np.zeros((N_ligs, 3, 3), dtype=float)
	Rz               = np.zeros((N_ligs, 3, 3), dtype=float)
	
	phi_random_idx = np.random.randint(angle_sampling, size=N_ligs)
	phi_random     = phi[phi_random_idx] * DEG2RAD
	
	Rx[:,0,0]        = 1.
	Rx[:,1,1]        =  np.cos(phi_random)
	Rx[:,2,2]        =  np.cos(phi_random)
	Rx[:,1,2]        = -np.sin(phi_random)
	Rx[:,2,1]        =  np.sin(phi_random)
	    
	phi_random_idx = np.random.randint(angle_sampling, size=N_ligs)
	phi_random     = phi[phi_random_idx] * DEG2RAD
	
	Ry[:,1,1]        = 1.
	Ry[:,0,0]        =  np.cos(phi_random)
	Ry[:,2,2]        =  np.cos(phi_random)
	Ry[:,2,0]        = -np.sin(phi_random)
	Ry[:,0,2]        =  np.sin(phi_random)
	
	phi_random_idx = np.random.randint(angle_sampling, size=N_ligs)
	phi_random     = phi[phi_random_idx] * DEG2RAD
	
	Rz[:,2,2]        = 1.
	Rz[:,0,0]        =  np.cos(phi_random)
	Rz[:,1,1]        =  np.cos(phi_random)
	Rz[:,0,1]        = -np.sin(phi_random)
	Rz[:,1,0]        =  np.sin(phi_random)
	
	xyz_ligand_new = np.einsum('aij,akj->aki', Rx, xyz_ligand_origin)
	xyz_ligand_new = np.einsum('aij,akj->aki', Ry, xyz_ligand_new)
	xyz_ligand_new = np.einsum('aij,akj->aki', Rz, xyz_ligand_new)
	
	xyz_ligand_new = xyz_sphere + xyz_ligand_new

	print "Saving structures..."

	for i in range(N_ligs):
	
		pt.write_trajectory(prefix+"_%d."%i+args.outformat, xyz_ligand_new, top=lig.top, overwrite=True, frame_indices=[i])

if __name__ == "__main__":

	args = parser.parse_args()

	main()