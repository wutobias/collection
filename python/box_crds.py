### Numerical libraries
import numpy as np

### MD specific libraries
import mdtraj as md

### Use field classes from PyGIST
from io_operations import field
from io_operations import write_files ### Mainly for debugging

### Other
from helpers import are_you_numpy

class box_crds(object):

	def __init__(self, topo_path, bins, strc_path_ref=None, pdb_path_ref=None):

		self.topo               = md.load_topology(topo_path)
		self.solvent_O_idxs     = self.topo.select("water and name O")
		self.bins               = bins
		self.grid_resolution    = np.array([1.,1.,1.])
		self.xx                 = np.zeros(3, dtype=float)
		self.yy                 = np.zeros(3, dtype=float)
		self.zz                 = np.zeros(3, dtype=float)
		self.center             = np.zeros(3, dtype=float)
		self.f2r                = np.eye(3,3)
		self.f                  = None
		self.fitting            = False
		self.ref                = None
		self.ref_sele           = None
		self.sele               = None

		if strc_path_ref != None and pdb_path_ref != None:
			self.ref      = md.load_frame(strc_path_ref, 0, top=pdb_path_ref)
			self.ref_sele = self.ref.topology.select("name CA or name N or name C")
			self.sele     = self.topo.select("name CA or name N or name C")
			self.fitting  = True


	def set_frame(self, frame, unit_cell, center):

		self.frame       = frame
		if self.fitting:
			self.frame.superpose(self.ref,
                                 frame=0,
                                 atom_indices=self.sele,
                                 ref_atom_indices=self.ref_sele,
                                 parallel=True)

		# We must assume that incoming unit_cell has normalized principal vectors
		self.f2r    = unit_cell * self.grid_resolution

		self.center = center

		self.set_field()


	def set_field(self):

		### field constructor: __init__(self, Bins, Frac2Real=None, Delta=None, Origin=None, Center=None)
		self.f = field(Bins=self.bins, Frac2Real=self.f2r, Center=self.center)

	def get_inside_crds_frac(self):

		solvent_frac = self.f.get_frac(self.frame.xyz[0][self.solvent_O_idxs]*10.)
		_inside_idxs = np.where( (solvent_frac[:,0] >= .0) * (solvent_frac[:,0] <= self.bins[0]) * \
		                         (solvent_frac[:,1] >= .0) * (solvent_frac[:,1] <= self.bins[1]) * \
		                         (solvent_frac[:,2] >= .0) * (solvent_frac[:,2] <= self.bins[2]) )

		return solvent_frac[_inside_idxs]
