import pytraj as pt
import numpy as np
from scipy import spatial, stats

class watcor(object):

  def __init__(self, Traj, Parm, HS, start=0, stop=-1, step=1):
  
    self.traj_path  = Traj
    self.parm_path  = Parm
    self.hs_path    = HS

    self.traj      = pt.iterload(self.traj_path, self.parm_path, frame_slice=(start, stop, step))
    self.hs        = pt.load(self.hs_path)

  
  def gauss(self, x2, width2):
  
    y = 1./np.sqrt(2. * np.pi * width2) * np.exp(- x2/(2. * width2))
    return y

  def cor_pop(self, sele=":WAT@O", radius=1.0):

    mask = self.traj.top.select(sele)
    cor  = np.zeros((self.traj.n_frames, self.hs[0].n_atoms), dtype=float)
    radius2 = radius**2
  
    index=0
    for frame in self.traj.iterframe(autoimage=True):
  
      if index % 1000 == 0:

        print "Frame %s..." %index

      frame_pop = np.zeros(self.hs[0].n_atoms, dtype=float)
      dists2    = spatial.distance.cdist(self.hs[0].coordinates, frame[mask], metric="sqeuclidean")
      valids    = np.where(dists2 < radius2)
      
      frame_pop[valids[0]] = self.gauss(dists2[valids], radius2)
    
      cor[index] = np.copy(frame_pop)
      index += 1

    ### Each row of x represents a variable, and each column a single observation of all those variables
    return np.corrcoef(cor.T)
    
  def frac_pop(self, sele=":WAT@O", radius=1.0):

    mask = self.traj.top.select(sele)

    frac = np.zeros((self.hs[0].n_atoms), self.hs[0].n_atoms)
  
    for index, frame in enumerate(self.traj):
  
      if index % 100 == 0:

        print "Frame %s..." %index

      frame_pop  = np.zeros(len(self.hs[0].n_atoms))

      frame_pop[np.where\
               (spatial.distance.cdist\
               (self.hs[0].coordinates, frame[mask], metric="euclidean") < radius)[0]] = 1

      frac[index] = np.copy(frame_pop)   

  def pymol_out(self, cor, name="session.py"):
  
    r = 0.15

    __header = str()
    __header += "from pymol import cmd\n"
    __header += "from pymol.cgo import *\n"
    __header += "cmd.load(\"%s\", \"HS\")\n" %self.hs_path
    __cor08_10   = "obj=[]\n"
    __cor06_08   = "obj=[]\n"
    __cor04_06   = "obj=[]\n"
    __cor02_04   = "obj=[]\n"
    __corm08_m10 = "obj=[]\n"
    __corm06_m08 = "obj=[]\n"
    __corm04_m06 = "obj=[]\n"
    __corm02_m04 = "obj=[]\n"
    
    corr_data = "#HS1 HS2 corr\n"
    
    for row_index, row in enumerate(cor):
    
      for column_index, column in enumerate(row):

        ### Matrix is symmetric
        if column_index >= row_index:
          break
      
        if -0.2 < column < 0.2:
          continue
          
        corr_data += "%d %d %f\n" %(row_index, column_index, column)

        c = str()
        c += "obj.extend([CYLINDER, "
        ### Coordinates first atom
        c += "%s, %s, %s," %(self.hs[0].coordinates[row_index][0], self.hs[0].coordinates[row_index][1], self.hs[0].coordinates[row_index][2])
        ### Coordinates second atom
        c += "%s, %s, %s," %(self.hs[0].coordinates[column_index][0],self.hs[0].coordinates[column_index][1], self.hs[0].coordinates[column_index][2])
        ### Cylinder radius
        c += "%s," %r

        if 0.8 < column <= 1.0:
        
          __cor08_10 += c
          ### Color code
          __cor08_10 += "0.9961, 0.1550, 0.1550, 0.9961, 0.1550, 0.1550, ])\n"

        if 0.6 < column <= 0.8:
        
          __cor06_08 += c
          __cor06_08 += "0.9961, 0.3638 , 0.3638, 0.9961, 0.3638, 0.3638, ])\n"

        if 0.4 < column <= 0.6:
        
          __cor04_06 += c
          __cor04_06 += "0.9961, 0.5725, 0.5725, 0.9961, 0.5725, 0.5725, ])\n"

        if 0.2 < column <= 0.4:
        
          __cor02_04 += c
          __cor02_04 += "0.9961, 0.7812, 0.7812, 0.9961, 0.7812, 0.7812, ])\n"

        if -0.8 >= column > -1.0:
        
          __corm08_m10 += c
          ### Color code
          __corm08_m10 += "0.1550, 0.1550, 0.9961, 0.1550, 0.1550, 0.9961, ])\n"

        if -0.6 >= column > -0.8:
        
          __corm06_m08 += c
          __corm06_m08 += "0.3638, 0.3638, 0.9961, 0.3638, 0.3638, 0.9961, ])\n"

        if -0.4 >= column > -0.6:
        
          __corm04_m06 += c
          __corm04_m06 += "0.5331, 0.5331, 0.9961, 0.5331, 0.5331, 0.9961, ])\n"

        if -0.2 >= column > -0.4:
        
          __corm02_m04 += c
          __corm02_m04 += "0.7937, 0.7937, 0.9961, 0.7937, 0.7937, 0.9961, ])\n"

    __cor08_10   += "cmd.load_cgo(obj, \"cor08_10  \")\n"
    __cor06_08   += "cmd.load_cgo(obj, \"cor06_08  \")\n"
    __cor04_06   += "cmd.load_cgo(obj, \"cor04_06  \")\n"
    __cor02_04   += "cmd.load_cgo(obj, \"cor02_04  \")\n"
    __corm08_m10 += "cmd.load_cgo(obj, \"cor-08_-10\")\n"
    __corm06_m08 += "cmd.load_cgo(obj, \"cor-06_-08\")\n"
    __corm04_m06 += "cmd.load_cgo(obj, \"cor-04_-06\")\n"
    __corm02_m04 += "cmd.load_cgo(obj, \"cor-02_-04\")\n"

    o = open(name, "w")

    o.write(__header)
    o.write(__cor08_10  )
    o.write(__cor06_08  )
    o.write(__cor04_06  )
    o.write(__cor02_04  )
    o.write(__corm08_m10)
    o.write(__corm06_m08)
    o.write(__corm04_m06)
    o.write(__corm02_m04)

    o.close()
    
    o = open(name.replace(".py", ".")+"corr.dat", "w")
    o.write(corr_data)
    o.close()
    