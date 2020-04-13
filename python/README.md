This is a collection of python scripts. Most of them are written in python
2.7, the newer ones are python >3.6 .

## Scripts
These are scripts which can be used as standalone.

#### FixHostFile.py
This script renumbers and reorders the atoms in a mol2/pdb file of a cyclodextrin molecule and writes it out in pdb format.

#### get_kb.py
Calculate KL divergence for two datasets

#### calc_total_charge_float.py
Calculate partial charge sum and print as floating point number from mol2 file

#### calc_total_charge.py
Calculate partial charge sum and print as (rounded) integer from mol2 file

#### normalize_B.py
Normalize B factors in pdb file according to neighboring atoms (within cutoff distance)

#### traj-k-means.py
K-Means clustering in cartesian or pca space

#### filter_clashes.py
filter clashes of structure and query structure

#### fixparm.py
fix corrupt amber parameter files

#### ITC_fit.py
Fit ITC data from multiple titrations in different buffers and perform enthalpy correction

#### matrix_align.py
rotate/translate structures using an object matrix calculated from a superposition of some 3rd structure

#### make_supercell.py
Generate supercells from crystal structure pdb file which contains appropriate unit cell 
information (needs xtalmd module)

#### watcor_analysis.py
Calculate correlation between Hydration sites (needs watcor module)

#### clust_wat_frac_bruteforce.py
Cluster water molecules with brute force
 
#### make_pca.py
Perform pca on MD trajectory, can also calculate Z-score of overlap with another trajectory

#### entropy_trans.py
Calculate entropy of water molecules around some solute residue.

#### match_charges.py
match the charges from mol2 file A to mol2 file B. The molecules must be identical, but can have different atom ordering, atom names, etc...

#### packmeup.py
Extract ligands from ensemble/trajectory and place it randomly around a protein.

#### strip_random_wat.py
Strip water molecules randomly from already build simulation box and rewrite.

#### get_cog.py
Calcualte center of geometry for given residue selection

#### esptotbl.py
Write esp file to table

#### esptopdb.py
Write esp file to pdb

#### clust_wat.py
Cluster water molecules using a pre-calculated (e.g. with cpptraj) occupancy grid

#### mapconv.py
Can perform various calculations on grids obtained from GIST analysis. Developed during my master thesis in 2015.

#### molecular_diameter.py
Calculate molecular diameter.

## Classes and functions
These files contain classes and functions that can be helpful. Some of the scripts above depend on these.

#### plot_hist.py
Contains some high-level functions for plotting 1d and 2d histograms.

#### io_operations.py
Different methods for reading and writing dx and pdb files.

#### helpers.py
Various helper functions.

#### spatial.py
Some helper functions for spatial transformations.

#### watcor
Module containing methods for calculating correlations between Hydration sites.

#### xtalmd
Unfinished module for preparing and running MD simulations of crystal lattices. Used by make_supercell.py

#### box_crds.py
Handling atom fractional coordinates in simulation box
