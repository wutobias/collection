import pytraj as pt
import numpy as np
from scipy import spatial, integrate
import matplotlib.pyplot as plt

class occ_analysis(object):

  def __init__(self, Traj, Parm, HS, start=0, stop=-1, step=1, mask=":WAT@O", radius=1.5, dt=1.0):

    __doc__="""

    Traj  : Path to trajectory file
    Parm  : Path to parameter file
    HS  : Path to cluster center (hydration sites) file
    start : Start frame trajectory, default 0
    stop  : Stop frame trajectory, default -1 (last)
    step  : Stride every step frame, default 1
    mask  : selection mask used for water calculation, default \":WAT@O\"
    radius: radius of cluster center
    dt  : integration step for autocorrelation in pico seconds
    """

    self.traj_path = Traj
    self.parm_path = Parm
    self.hs_path   = HS
    self.mask      = mask
    self.radius    = radius
    self.dt        = dt

    self.traj    = pt.iterload(self.traj_path, self.parm_path, frame_slice=(start, stop, step), mask=self.mask)
    self.hs      = pt.load(self.hs_path)

    ### This contains all data of the time series
    ### rows contain frames, columns contain different hydration sites
    self.series    = np.zeros((self.traj.n_frames, self.hs.n_atoms), dtype=np.int_)
    self.autocorr  = np.zeros((self.traj.n_frames, self.hs.n_atoms))
    self.crosscorr = np.zeros((self.traj.n_frames, self.hs.n_atoms, self.hs.n_atoms))

    self.__get_data()


  def __get_data(self):

    mask = self.traj.top.select(self.mask)

    for index, frame in enumerate(self.traj):
  
      if index % 1000 == 0:

        print "Frame %s..." %index

      ### Store indices of water molecules
      ### in hydration sites in current frame. 
      __frame_idx     = -1 * np.ones(self.hs.n_atoms)
      __hs_idx, __wat_idx = np.where\
                  (spatial.distance.cdist\
                  (self.hs[0].coordinates, frame[mask], metric="euclidean") < self.radius)

      __frame_idx[__hs_idx] = mask[__wat_idx]

      self.series[index]  = np.copy(__frame_idx)

    ### Store whether HS Y contains water at frame X (=1)
    ### or not (=0).
    self.__binary_occupied = np.zeros((self.traj.n_frames, self.hs.n_atoms))
    self.__binary_occupied[np.where(self.series>=0)] = 1.0


  def __pettitt(self,k):

    ### Autocorrelation functional according to Pettitt.
    ### See DOI: 10.1016/S0006-3495(00)76533-7

    ### Integration sample points t:
    sample_t            = np.arange(self.traj.n_frames-k)

    autocor             = np.zeros((self.traj.n_frames-k, self.hs.n_atoms))
    __value_t, __hs     = np.where(self.series[sample_t]==self.series[sample_t+k])
    __hs_self           = np.where(__hs==__hs)
    autocor[__value_t]  = self.__binary_occupied[__value_t, __hs_self] * self.__binary_occupied[__value_t+k, __hs_self]
    autocor             = np.sum(autocor, axis=0)

    ### Normalize autocorrelation
    autocor            /= np.sum(self.__binary_occupied, axis=0)

    return autocor


  def autocorrelate(self):

    for k in range(self.traj.n_frames):

      self.autocorr[k] = self.__pettitt(k)


  def print_ac(self, HS_idx=None, prefix="autocor"):

    if HS_idx == None:

      HS_idx = range(self.hs.n_atoms)

    for hs in HS_idx:

      plt.plot(np.arange(self.traj.n_frames,dtype=float)*self.dt, self.autocorr[:,hs])
      axes = plt.gca()
      axes.set_xlim([0,self.traj.n_frames])
      axes.set_ylim([0,1])
      plt.ylabel("R(t)")
      plt.xlabel("t[ps]")
      plt.title("Hydration Site %d" %hs)
      plt.savefig("HS_%d_autocor.png" %hs)
      plt.close('all')

  def __my_cc(self, k):

    ### Integration sample points t:
    sample_t            = np.arange(self.traj.n_frames-k)

    crosscor            = np.zeros((self.traj.n_frames-k, self.hs.n_atoms, self.hs.n_atoms))
    __value_t, __hs     = np.where(self.series[sample_t]==self.series[sample_t+k])

    for idx1 in range(self.hs.n_atoms):

      for idx2 in range(self.hs.n_atoms):

        if idx1 < idx2:

          __hs_1        = np.where(__hs==id1)
          __hs_2        = np.where(__hs==id2)
          crosscor[__value_t]  = self.__binary_occupied[__value_t, __hs_1] * self.__binary_occupied[__value_t+k, __hs_2]
          crosscor             = np.sum(crosscor, axis=0)

    ### Normalize autocorrelation
    autocor            /= np.sum(self.__binary_occupied, axis=0)

    crosscor = np.zeros((self.traj.n_frames-k, self.hs.n_atoms, self.hs.n_atoms))

    for idx1 in range(self.hs.n_atoms):

      for idx2 in range(self.hs.n_atoms):

        if idx1 < idx2:

          __value_t, value_tau = np.where(self.series[sample_t]==self.series[sample_t+k,idx2])

          crosscor[__value_t, idx1, idx2] = self.__binary_occupied[__value_t, idx1] * self.__binary_occupied[value_tau, idx2]
        
    crosscor = np.sum(crosscor, axis=0)

    ### Normalize
    for idx1 in range(self.hs.n_atoms):

      for idx2 in range(self.hs.n_atoms):

        if idx1 < idx2:

          crosscor[idx1,idx2] /= np.sum(self.__binary_occupied[:,idx1]*self.__binary_occupied[:,idx2], axis=0)

    return crosscor


  def cross_correlate(self):

    for k in range(self.traj.n_frames):

      self.crosscorr[k] = self.__my_cc(k)

