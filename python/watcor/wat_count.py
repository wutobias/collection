import pytraj as pt
import numpy as np
from scipy import spatial, stats

class wat_count(object):

  def __init__(self, Traj, Parm, HS, start=0, stop=-1, step=1, mask=":WAT@O", radius=1.5):

    __doc__="""

    Traj  : Path to trajectory file
    Parm  : Path to parameter file
    HS    : Path to cluster center (hydration sites) file
    start : Start frame trajectory, default 0
    stop  : Stop frame trajectory, default -1 (last)
    step  : Stride every step frame, default 1
    mask  : selection mask used for water calculation, default \":WAT@O\"
    radius: radius of cluster center
    """

    self.traj_path = Traj
    self.parm_path = Parm
    self.hs_path   = HS
    self.mask      = mask
    self.radius    = radius

    self.traj      = pt.iterload(self.traj_path, self.parm_path, frame_slice=(start, stop, step), mask=self.mask)
    self.hs        = pt.load(self.hs_path)

    ### This contains all data of the time series
    ### rows contain frames, columns contain different hydration sites
    self.count    = np.zeros(self.hs.n_atoms)
    self.occ      = np.zeros(self.hs.n_atoms)

    self.__get_data()

  def __get_data(self):

    mask = self.traj.top.select(self.mask)

    for index, frame in enumerate(self.traj):
  
      if index % 1000 == 0:

        print "Frame %s..." %index

      ### Store indices of water molecules
      ### in hydration sites in current frame. 
      self.count[np.where\
                  (spatial.distance.cdist\
                  (self.hs[0].coordinates, frame[mask], metric="euclidean") < self.radius)[0]] += 1

    self.occ = self.count / self.traj.n_frames

  def print_count(self):

    for hs_idx in range(self.hs.n_atoms):

      print "HS %d: Count %d Occ %6.3f" %(hs_idx, self.count[hs_idx], self.occ[hs_idx])
