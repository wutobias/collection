#!/usr/bin/env python

# In[132]:

import numpy as np
import pytraj as pt
import scipy.spatial as scipy_spatial
import sys
from collections import OrderedDict as dict

DEG2RAD        = np.pi/180.
RAD2DEG        = 180./np.pi
N_A            = 6.02214 * 10**23
### R_gas is unit kcal*K^-1*mol^-1*
R_gas          = 1.9872036*10**-3
TEMP           = 300.

if len(sys.argv) != 4:

    print "Usage: entropy_AR.py [prmtop] [traj] [solute_resid]"
    exit()

start = 1.0
stop  = 5.0
step  = 0.2

parm_path = sys.argv[1]
traj_path = sys.argv[2]
sol_mask  = ":%s&!@H=" %(sys.argv[3])


# In[50]:

class entropy(object):
    
    def __init__(self, X):

        self.X      = X
        self.N      = X.shape[0]
        
        self.quart  = False
        
        if self.X.shape[1] == 4:
            
            self.quart = True

        self.kdtree = scipy_spatial.cKDTree(self.X)
        
    def get_geodesic(self):
        
        dd = np.zeros(self.N, dtype=float)
        ii = np.zeros(self.N, dtype=int)
        
        for i in range(self.N):
        
            q_dot       = np.einsum('ij,j->i', self.X, self.X[i])
            q_dot_upper = np.where(q_dot> 1.0)
            q_dot_lower = np.where(q_dot<-1.0)
        
            q_dot[q_dot_upper] =  1.0
            q_dot[q_dot_lower] = -1.0
            
            sele       = np.arange(self.N-1)
            sele[i:]  += 1
            q_dist     = 2.*np.arccos(q_dot)[sele]
            
            _min       = np.argmin(q_dist)
            
            dd[i]      = q_dist[_min]
            ii[i]      = sele[_min]
            
        return dd,ii
    
    def get_normofdifference(self):
        
        dd = np.zeros(self.N, dtype=float)
        ii = np.zeros(self.N, dtype=int)
        
        for i in range(self.N):
        
            sele        = np.arange(self.N-1)
            sele[i:]   += 1

            q_norm      = np.zeros(self.N-1, dtype=float)
            
            q_diff_1    = self.X[i] - self.X
            q_diff_2    = self.X[i] + self.X
            q_norm_1    = np.linalg.norm(q_diff_1, axis=1)[sele]
            q_norm_2    = np.linalg.norm(q_diff_2, axis=1)[sele]
            
            q_norm_min                 = np.where(q_norm_1>q_norm_2)
            q_norm_min_inv             = np.ones(self.N-1, dtype=bool)
            q_norm_min_inv[q_norm_min] = False
            q_norm[q_norm_min]         = q_norm_2[q_norm_min]
            q_norm[q_norm_min_inv]     = q_norm_1[q_norm_min_inv]

            _min                       = np.argmin(q_norm)

            dd[i]                      = 2.*q_norm[_min]
            ii[i]                      = sele[_min]
            
        return dd,ii        
        
    
    def get_knn_eucl(self,k=1):
        
        dd, ii = self.kdtree.query(self.X, k+1)
        
        return dd[:,1:],ii[:,1:]
    
    def get_knn_vol(self,k=1):
        
        if self.quart:
            
            #dd, ii = self.get_geodesic()
            dd, ii = self.get_normofdifference()
        
        else:
        
            dd, ii = self.get_knn_eucl(k)
        
        return 4./3.*np.pi*np.power(dd, 3.)
    
    def get_pdf(self, k=1):
        
        v = self.get_knn_vol(k)
        
        return 8.*np.pi**2/v
    
    def get_entropy(self, k=1):
        
        pdf = self.get_pdf(k=1)
        
        return R_gas/self.N * (np.sum(np.log(pdf)) + np.euler_gamma)


# In[133]:

traj       = pt.iterload(traj_path, parm_path)
O_idx      = pt.select_atoms(traj.top, ':WAT@O')
Solute_idx = pt.select_atoms(traj.top, sol_mask)
N_frames   = traj.n_frames

### Key : distance
### Item: [water_O_idx, frame_idx, count]

crd_dict   = dict()

for d in range(int((stop-start)/step)):
    
    crd_dict[start+d*step] = list()

for frame in traj.iterframe():

    pt.autoimage(frame, top=traj.topology)
    solute  = frame[Solute_idx]
    solvent = frame[O_idx]
    
    dists   = scipy_spatial.distance.cdist(solvent, solute)
    for dist, l in crd_dict.items():
        
        valids = np.where(dists < dist)[0]
        valids = np.unique(valids)
        
        N_valids = valids.shape[0]
        
        if N_valids > 0:
        
            for N_i in range(N_valids):
            
                l.append(solvent[valids[N_i]])
                
    del dists

print "# dist TdS[kcal/mol]"
    
for dist, l in crd_dict.items():
    
    L = len(l)
    
    if L>0:
    
        crds = np.array(l)
        e = entropy(crds)
        print dist, TEMP*e.get_entropy()
    else:
        
        print dist, "0"

