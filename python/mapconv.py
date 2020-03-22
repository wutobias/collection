#!/usr/bin/env python

'''
#################################################################
#                                                               #
#  mapconv.py is written by Tobias Wulsdorf and comes           #
#  with no warrenty. This python-script can do different        #
#  operations on density maps obtained from molecular dynamics  #
#  simulations. The main focus however are fluid density and    #
#  interaction maps.                                            #
#  Type "mapconv.py --help" for a short overview of all         #
#  functionalities.                                             #
#                                                               #
#     e-mail: tobias.wulsdorf@pharmazie.uni-marburg.de          #
#                                                               #
# System requirements:                                          #
#  python-version > 2.7.3                                       #
#                                                               #
# Python requirements:                                          #
#  os                                                           #
#  sys                                                          #
#  math                                                         #
#  argparse                                                     #
#  time                                                         #
#  MDAnalysis                                                   #
#  subprocess                                                   #
#  numpy 1.9.2                                                  #
#  scipy                                                        #
#  random                                                       #
#  Bio.PDB                                                      #
#  collections                                                  #
#  tempfile                                                     #
#                                                               #
# Others:                                                       #
#  msms                                                         #
#                                                               #
#################################################################
'''

'''
## G e n e r a l   i n f o r m a t i o n ##

- All output obtained from Grid Inhomogenous Solvation Theory (GIST) calculations
  is given in terms of a displacement process:

  water_rec --> water_bulk

- GIST displacement thermodynamics are given as volume normalized (NOT density normalized!)

'''

import os
import sys
import math
import argparse
import time
import MDAnalysis
import subprocess
import numpy as np
import scipy.optimize
from scipy import ndimage
from scipy.spatial import *
import scipy.stats
import scipy.interpolate
import scipy.special
import random
from collections import defaultdict
from collections import OrderedDict
import tempfile
from Bio.PDB import Selection
from Bio.PDB.AbstractPropertyMap import AbstractPropertyMap
from Bio.PDB.Polypeptide import is_aa

path2script = sys.argv[0]

parser = argparse.ArgumentParser(prog='mapconv')


#Constants
GAS_CONST_KCAL = 0.0019872041

###Functions###

def read_map (path, format, read_header=False, gist_data_add={}, reference_density_add=None):

   # {[x_coord,y_coord,z_coord]:value}
   dict_frac = {}

   origin_x = None
   n_x = None
   min_x = None
   max_x = None
   dim_x = None
   alpha = None
   bins_x = None

   origin_y = None
   n_y = None
   min_y = None
   max_y = None
   dim_y = None
   beta = None
   bins_y = None

   origin_z = None
   n_z = None
   min_z = None
   max_z = None
   dim_z = None
   gamma = None
   bins_z = None

   map_file_ref = open(path,"r")
   map_file = map_file_ref.readlines()
   map_file_ref.close()

   if (format == 'xplor'):

      ### FOR XPLOR FORMAT ####
   
      # read header information
      
      start_row = -1

      for i, item in enumerate(map_file):
      
        if len(item.rstrip().split()) > 0 and ( item.rstrip().split()[0] == 'ZYX' or  item.rstrip().split()[0] == 'XYZ'):
        
          start_row = i - 2
        
          break
   
      n_x = float(map_file[start_row].rstrip().split()[0])
      min_x = int(map_file[start_row].rstrip().split()[1])
      max_x = int(map_file[start_row].rstrip().split()[2])
      dim_x = float(map_file[1+start_row].rstrip().split()[0])
      alpha = float(map_file[1+start_row].rstrip().split()[3])
      bins_x = max_x - min_x + 1
   
      n_y = float(map_file[start_row].rstrip().split()[3])
      min_y = int(map_file[start_row].rstrip().split()[4])
      max_y = int(map_file[start_row].rstrip().split()[5])
      dim_y = float(map_file[1+start_row].rstrip().split()[1])
      beta = float(map_file[1+start_row].rstrip().split()[4])
      bins_y = max_y - min_y + 1
   
      n_z = float(map_file[start_row].rstrip().split()[6])
      min_z = int(map_file[start_row].rstrip().split()[7])
      max_z = int(map_file[start_row].rstrip().split()[8])
      dim_z = float(map_file[1+start_row].rstrip().split()[2])
      gamma = float(map_file[1+start_row].rstrip().split()[5])
      bins_z = max_z - min_z + 1

      cos_a = ( math.cos(math.pi*beta/180) * math.cos(math.pi*gamma/180) - math.cos(math.pi*alpha/180) ) / ( math.sin(math.pi*beta/180) *  math.sin(math.pi*gamma/180))
      sin_a = math.sqrt(1.0 - cos_a**2)

      origin_x = dim_x  / n_x * min_x  +  dim_y / n_y * math.cos(math.pi*gamma/180) * min_y  +  dim_z / n_z * math.cos(math.pi*beta/180) * min_z
      origin_y = 0                     +  dim_y / n_y * math.sin(math.pi*gamma/180) * min_y  -  dim_z / n_z * math.sin(math.pi*beta/180) * cos_a * min_z
      origin_z = 0                     +  0                                                  +  dim_z / n_z * math.sin(math.pi*beta/180) * sin_a * min_z

      origin = []
   
      origin.append(origin_x)
      origin.append(origin_y)
      origin.append(origin_z)
         
      n = []
      
      n.append(n_x)
      n.append(n_y)
      n.append(n_z)
         
      dim = []
   
      dim.append(dim_x)
      dim.append(dim_y)
      dim.append(dim_z)
   
      bins = []
   
      bins.append(bins_x)
      bins.append(bins_y)
      bins.append(bins_z)

      if read_header:

        dict_map = {'map':dict_frac, 'origin':origin, 'n':n, 'dim':dim, 'bins':bins, 'alpha':alpha, 'beta':beta, 'gamma':gamma, 'format':format}

        return dict_map
   
      sector_factor = bins_x * bins_y / 6
      
      if ((bins_x * bins_y) % 6 != 0):
         
         sector_factor = sector_factor + 1
                          
      # read in map

      for z_index, z_coord in enumerate(range(0, bins_z)):
   
        sector = []
      
        for line in map_file[4+start_row + z_index * sector_factor + z_index : 4+start_row + (z_index + 1) * sector_factor + z_index]: 

          for line_element in range(len(line.rstrip()) / 12):

            sector.append(line.rstrip()[line_element*12 : (line_element + 1) * 12])

        for y_index, y_coord in enumerate(range(0, bins_y)):

          for x_index, x_coord in enumerate(range(0, bins_x)):

            dict_frac[x_coord,y_coord,z_coord] = float(sector[x_index + y_index * bins_x])


      print "\nFinished reading xplor-map."

      dict_map = {'map':dict_frac, 'origin':origin, 'n':n, 'dim':dim, 'bins':bins, 'alpha':alpha, 'beta':beta, 'gamma':gamma, 'format':format}
   
      return dict_map      


   if (format == 'dx'):

      start_row = -1

      for i, item in enumerate(map_file):

        if len(item.rstrip().split()) > 1 and item.rstrip().split()[0] == 'object' and item.rstrip().split()[1] == '1':

          start_row = i

      origin_x = float(map_file[start_row + 1].rstrip().split()[1])
      bins_x = int(map_file[start_row].rstrip().split()[5])
      n_x = bins_x 
      dim_x = bins_x * float(map_file[start_row+2].rstrip().split()[1])
      alpha = 90.0
   
      origin_y = float(map_file[start_row + 1].rstrip().split()[2])
      bins_y = int(map_file[start_row].rstrip().split()[6])
      n_y = bins_y 
      dim_y = bins_y * float(map_file[start_row+3].rstrip().split()[2])
      beta = 90.0

      origin_z = float(map_file[start_row + 1].rstrip().split()[3])
      bins_z = int(map_file[start_row].rstrip().split()[7])
      n_z = bins_z
      dim_z = bins_z * float(map_file[start_row+4].rstrip().split()[3])
      gamma = 90.0

      origin = []
   
      origin.append(origin_x)
      origin.append(origin_y)
      origin.append(origin_z)
         
      n = []
      
      n.append(n_x)
      n.append(n_y)
      n.append(n_z)
         
      dim = []
   
      dim.append(dim_x)
      dim.append(dim_y)
      dim.append(dim_z)
   
      bins = []
   
      bins.append(bins_x)
      bins.append(bins_y)
      bins.append(bins_z)

      if read_header:

        dict_map = {'map':dict_frac, 'origin':origin, 'n':n, 'dim':dim, 'bins':bins, 'alpha':alpha, 'beta':beta, 'gamma':gamma, 'format':format}

        return dict_map

      dx_list = []

      rows_of_data = (bins_x * bins_y * bins_z) / 3

      if (bins_x * bins_y * bins_z) % 3 != 0:

        rows_of_data = rows_of_data + 1

      for i in range(start_row+7, start_row+7+rows_of_data):
         for j in map_file[i].rstrip().split():

            dx_list.append(j)

      for x in range(0, bins_x):
         for y in range(0, bins_y):
            for z in range(0, bins_z):

               dict_frac[x,y,z] = float(dx_list[x * bins_y * bins_z  + y * bins_z + z])
               
      print "\nFinished reading dx-map."

      dict_map = {'map':dict_frac, 'origin':origin, 'n':n, 'dim':dim, 'bins':bins, 'alpha':alpha, 'beta':beta, 'gamma':gamma, 'format':format}
   
      return dict_map      



   if (format == 'Eww_ij'):

    addbins = gist_data_add['bins']
    gO_add  = gist_data_add['gO']
    dim_add = gist_data_add['dim']
    n_add   = gist_data_add['n']
    dV_add  = (dim_add[0] / n_add[0]) * (dim_add[1] / n_add[1]) * (dim_add[2] / n_add[2])

    Eww_ij    = defaultdict(dict)

    for line in map_file:

      i   = int(line.rstrip().split()[0])
      j   = int(line.rstrip().split()[1])

      x_i =  i / (addbins[1] * addbins[2])
      y_i = (i % (addbins[1] * addbins[2])) / addbins[2]
      z_i = (i % (addbins[1] * addbins[2])) % addbins[2]

      x_j =  j / (addbins[1] * addbins[2])
      y_j = (j % (addbins[1] * addbins[2])) / addbins[2]
      z_j = (j % (addbins[1] * addbins[2])) % addbins[2]

#      Eww_ij[x_i, y_i, z_i][x_j, y_j, z_j] = float(line.rstrip().split()[2]) / (gO_add[x_j, y_j, z_j]*dV_add*reference_density_add)
#      Eww_ij[x_j, y_j, z_j][x_i, y_i, z_i] = float(line.rstrip().split()[2]) / (gO_add[x_i, y_i, z_i]*dV_add*reference_density_add)

      Eww_ij[x_i, y_i, z_i][x_j, y_j, z_j] = float(line.rstrip().split()[2])
      Eww_ij[x_j, y_j, z_j][x_i, y_i, z_i] = float(line.rstrip().split()[2])

    print "Finished reading Eww_ij matrix. Length %d" %len(Eww_ij)

    return Eww_ij

   if (format == 'gist_out'):

      start_row = -1

      for i, item in enumerate(map_file):

        if len(item.rstrip().split()) > 1 and item.rstrip().split()[0] == 'voxel' and item.rstrip().split()[1] == 'xcoord':

          start_row = i + 1

          break

      bins_z = -1
      bins_y = -1
      bins_x = -1

      dim_z = 0.0
      dim_y = 0.0
      dim_x = 0.0

#in order to get the information from the gist-out file into a format that can be handled by this python-script, one needs to find out the bin-dimensions first...

      z_start = map_file[start_row].rstrip().split()[3]
      y_start = map_file[start_row].rstrip().split()[2]
      x_start = map_file[start_row].rstrip().split()[1]

      found_bins_z = False
      found_bins_y = False
      found_bins_x = False

      for i, line in enumerate(map_file[start_row+1:]):

        if found_bins_z and not found_bins_y and line.rstrip().split()[2] == y_start:

          bins_y = (i+1) / bins_z

          break

        if not found_bins_z and line.rstrip().split()[3] == z_start:

          bins_z = i + 1

          found_bins_z = True

      bins_x = (len(map_file) - start_row) / (bins_z * bins_y)

      n_x = bins_x
      n_y = bins_y
      n_z = bins_z

      dim_z = abs(float(z_start) - float(map_file[start_row + 1].rstrip().split()[3])) * float(bins_z)
      dim_y = abs(float(y_start) - float(map_file[start_row + 1 + bins_z].rstrip().split()[2])) * float(bins_y)
      dim_x = abs(float(x_start) - float(map_file[start_row + 1 + bins_y * bins_z].rstrip().split()[1])) * float(bins_x)

      origin_z = float(map_file[start_row].rstrip().split()[3]) - dim_z/float(bins_z) * 0.5
      origin_y = float(map_file[start_row].rstrip().split()[2]) - dim_y/float(bins_y) * 0.5
      origin_x = float(map_file[start_row].rstrip().split()[1]) - dim_x/float(bins_x) * 0.5

      alpha = 90.0
      beta  = 90.0
      gamma = 90.0

      origin = []

      origin.append(origin_x)
      origin.append(origin_y)
      origin.append(origin_z)

      n = []

      n.append(n_x)
      n.append(n_y)
      n.append(n_z)

      dim = []

      dim.append(dim_x)
      dim.append(dim_y)
      dim.append(dim_z)

      bins = []

      bins.append(bins_x)
      bins.append(bins_y)
      bins.append(bins_z)

#now we read all the gist-specific data...

      crd ={}
      pop = {}
      gO = {}
      gH = {}
      dTStrans_dens = {}
      dTStrans_norm = {}
      dTSorient_dens = {}
      dTSorient_norm = {}
      Esw_dens = {}
      Esw_norm = {}
      Eww_dens = {}
      Eww_norm_unref = {}
      Dipole_x_dens = {}
      Dipole_y_dens = {}
      Dipole_z_dens = {}
      Dipole_dens = {}
      neighbor_dens = {}
      neighbor_norm = {}
      order_norm = {}

      cos_a = ( math.cos(math.pi*beta/180) * math.cos(math.pi*gamma/180) - math.cos(math.pi*alpha/180) ) / ( math.sin(math.pi*beta/180) *  math.sin(math.pi*gamma/180))
      sin_a = math.sqrt(1.0 - cos_a**2)

      for x in range(0,bins_x):
        for y in range(0,bins_y):
          for z in range(0,bins_z):

            #Those coordinates are shifted.....narf!!!!!
            #crd[x,y,z]            = [float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[1]), float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[2]), float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[3])]

            X_coord               = origin[0] + dim[0] / n[0] * x + dim[1] / n[1] * math.cos(math.pi*gamma/180) * y + dim[2] / n[2] * math.cos(math.pi*beta/180) * z
            Y_coord               = origin[1] + 0                 + dim[1] / n[1] * math.sin(math.pi*gamma/180) * y - dim[2] / n[2] * math.sin(math.pi*beta/180) * cos_a * z
            Z_coord               = origin[2] + 0                 + 0                                               + dim[2] / n[2] * math.sin(math.pi*beta/180) * sin_a * z
            crd[x,y,z]            = [X_coord, Y_coord, Z_coord]

            pop[x,y,z]            = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[4])
            gO[x,y,z]             = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[5])
            gH[x,y,z]             = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[6])
            dTStrans_dens[x,y,z]  = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[7])
            dTStrans_norm[x,y,z]  = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[8])
            dTSorient_dens[x,y,z] = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[9])
            dTSorient_norm[x,y,z] = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[10])
            Esw_dens[x,y,z]       = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[11])
            Esw_norm[x,y,z]       = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[12])
            Eww_dens[x,y,z]       = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[13])
            Eww_norm_unref[x,y,z] = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[14])
            Dipole_x_dens[x,y,z]  = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[15])
            Dipole_y_dens[x,y,z]  = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[16])
            Dipole_z_dens[x,y,z]  = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[17])
            Dipole_dens[x,y,z]    = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[18])
            neighbor_dens[x,y,z]  = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[19])
            neighbor_norm[x,y,z]  = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[20])
            order_norm[x,y,z]     = float(map_file[start_row + x * bins_y * bins_z  + y * bins_z + z].rstrip().split()[21])

      dict_map = {'map':gO, 'origin':origin, 'n':n, 'dim':dim, 'bins':bins, 'alpha':alpha, 'beta':beta, 'gamma':gamma, 'format':format, 'crd':crd, 'pop':pop, 'gO':gO, 'gH':gH, \
      'dTStrans_dens':dTStrans_dens, 'dTStrans_norm':dTStrans_norm, 'dTSorient_dens':dTSorient_dens, 'dTSorient_norm': dTSorient_norm, 'Esw_dens':Esw_dens, 'Esw_norm':Esw_norm, \
      'Eww_dens': Eww_dens, 'Eww_norm_unref':Eww_norm_unref, 'Dipole_x_dens':Dipole_x_dens, 'Dipole_y_dens':Dipole_y_dens, 'Dipole_z_dens':Dipole_z_dens, 'Dipole_dens': Dipole_dens, \
      'neighbor_dens':neighbor_dens, 'neighbor_norm':neighbor_norm, 'order_norm':order_norm}

      print "Finished reading gist-out file."

      return dict_map 

def crop_Eww_matrix(path, gist_data_add, crop_grid=None):
#If crop_grid is empty we read the whole Eww_matrix

    if crop_grid==None:
      crop_grid=gist_data_add['crd']

    addbins = gist_data_add['bins']
    gO_add  = gist_data_add['gO']
    dim_add = gist_data_add['dim']
    n_add   = gist_data_add['n']
    dV_add  = (dim_add[0] / n_add[0]) * (dim_add[1] / n_add[1]) * (dim_add[2] / n_add[2])

    Eww  = 0.0

    map_file_temp = open(path, "r")
    map_file      = map_file_temp.readlines()
    map_file_temp.close()
    #with open(path) as map_file:
    for line in map_file:

      i   = int(line.rstrip().split()[0])
      j   = int(line.rstrip().split()[1])

      x_i =  i / (addbins[1] * addbins[2])
      y_i = (i % (addbins[1] * addbins[2])) / addbins[2]
      z_i = (i % (addbins[1] * addbins[2])) % addbins[2]

      x_j =  j / (addbins[1] * addbins[2])
      y_j = (j % (addbins[1] * addbins[2])) / addbins[2]
      z_j = (j % (addbins[1] * addbins[2])) % addbins[2]

      if [x_i, y_i, z_i] in crop_grid.keys() or [x_j, y_j, z_j] in crop_grid.keys():
        Eww += float(line.rstrip().split()[2])

    return Eww


#The following is copied from Bio.PDB package.
#The alteration of atomic radii in get_MS is new

def _read_vertex_array(filename):
    """
    Read the vertex list into a Numeric array.
    """
    with open(filename, "r") as fp:
        vertex_list=[]
        for l in fp.readlines():
            sl=l.split()
            if not len(sl)==9:
                # skip header 
                continue
            vl = [float(x) for x in sl[0:3]]
            vertex_list.append(vl)
    return np.array(vertex_list)


def get_MS(pdb_file, radius_add=1.0, PDB_TO_XYZR="pdb_to_xyzr", MSMS="msms"):
    """
    Return a Numeric array that represents
    the vertex list of the molecular surface.

    PDB_TO_XYZR --- pdb_to_xyzr executable (arg. to os.system)
    MSMS --- msms executable (arg. to os.system)
    """
    # extract xyz and set radii
    xyz_tmp=tempfile.mktemp()
    PDB_TO_XYZR=PDB_TO_XYZR+" %s > %s"
    make_xyz=PDB_TO_XYZR % (pdb_file, xyz_tmp)
    os.system(make_xyz)
    assert os.path.isfile(xyz_tmp), \
        "Failed to generate XYZR file using command:\n%s" % make_xyz
    #change radii by radius_add
    xyz_altered=tempfile.mktemp()
    alter_by_radius_add = 'cat %s | awk \'BEGIN {OFS=\"\\t\"} {$4=$4+%f ; print}\' >  %s' %(xyz_tmp, radius_add, xyz_altered)
    os.system(alter_by_radius_add)
    # make surface
    surface_tmp=tempfile.mktemp()
    MSMS=MSMS+" -probe_radius 1.5 -if %s -of %s > "+tempfile.mktemp()
    make_surface=MSMS % (xyz_altered, surface_tmp)
    os.system(make_surface)
    surface_file=surface_tmp+".vert"
    assert os.path.isfile(surface_file), \
        "Failed to generate surface file using command:\n%s" % make_surface
    # read surface vertices from vertex file
    surface=_read_vertex_array(surface_file)
    # clean up tmp files
    # ...this is dangerous
    # os.system("rm "+xyz_tmp)
    # os.system("rm "+surface_tmp+".vert")
    # os.system("rm "+surface_tmp+".face")
    return surface


def gist_sphere_expander(water_positions, gist_data, Eww_matrix, reference_energy, calculation_density, reference_density, output, radius_add_end, radius_dr):

  radius  = 0.25
  summary = OrderedDict()
  output_header  = 'radius  '
  while radius <= radius_add_end:
    summary[radius] = gist_water_analysis (water_positions, gist_data, radius, Eww_matrix, reference_energy, calculation_density, reference_density, output, return_it=True, radius_start=0.0)
#write the header for output
    output_header += '%3.2f ' %radius
    radius += radius_dr

  output_header += '\n'
  key_list = summary[0.5].keys()
  
  for index, water_position in enumerate(water_positions.keys()):

    output_data = ''

    for key in key_list:
      output_data += '%s ' %key
      for radius in summary.keys():
        output_data += '%6.3f ' %summary[radius][key][water_position]
      output_data += '\n'
#write the file
    o = open('%s_%d_gist_sphere_analysis_expander.dat' %(output,index), "w")
    o.write(output_header+output_data)
    o.close()
  print ("Done")

def gist_cavity_expander(path2molecule, gist_data, Eww_matrix, reference_energy,  calculation_density, reference_density, output, radius_add_end, radius_dr, happy_only=False, unhappy_only=False):

  radius  = 0.0
  summary = OrderedDict()
  output_pr  = 'radius   '
  while radius <= radius_add_end:
    summary[radius] = gist_cavity(path2molecule, gist_data, Eww_matrix, reference_energy,  calculation_density, reference_density, output, radius, happy_only, unhappy_only, return_values=True)
#write the header for output
    output_pr += '%3.2f ' %radius
    radius += radius_dr

  output_pr += '\n'
  key_list = summary[0.0].keys()
  for key in key_list:
    output_pr += '%s ' %key
    for radius in summary.keys():
      output_pr += '%6.3f ' %summary[radius][key]
    output_pr += '\n'

#write the file
  o = open('%s_gist_cavity_analysis_expander.dat' %output, "w")
  o.write(output_pr)
  o.close()
  print ("Done")

def gist_cavity(path2molecule, gist_data, Eww_matrix, reference_energy,  calculation_density, reference_density, output, radius_add, happy_only=False, unhappy_only=False, return_values=False):

  if happy_only and unhappy_only:
    print "Calculation of unhappy and happy hydration at the same time is not possible."
    exit()

  if happy_only:
    print "This is happy water only."

  if unhappy_only:
    print "This is unhappy water only."

  crd            = gist_data['crd']
  #dTStrans_norm  = gist_data['dTStrans_norm']
  #dTSorient_norm = gist_data['dTSorient_norm']
  #Esw_norm       = gist_data['Esw_norm']
  #Eww_norm       = gist_data['Eww_norm_unref']
  gO_GIST        = gist_data['gO']

  bins           = np.array(gist_data['bins'])
  n              = np.array(gist_data['n'])
  dim            = np.array(gist_data['dim'])
  voxel_volume   = (dim[0] / n[0]) * (dim[1] / n[1]) * (dim[2] / n[2])

#correct gO and calculate Entropy_trans grid
  gO_sw          = {}
  dTStrans_norm  = {}
  dTSorient_norm = {}
  Esw_norm       = {}
  Eww_norm       = {}
  n_k            = {}

  for gO_k, gO_k_value in gO_GIST.items():
    
      gO_sw[gO_k]          = gO_k_value  * calculation_density / reference_density
      n_k  [gO_k]          = gO_sw[gO_k] * reference_density * voxel_volume
      
      dTSorient_norm[gO_k] = gist_data['dTSorient_norm'][gO_k] * n_k[gO_k] / voxel_volume
      Esw_norm      [gO_k] = gist_data['Esw_norm'][gO_k] * n_k[gO_k] / voxel_volume
      Eww_norm      [gO_k] = (reference_energy - 2.0*gist_data['Eww_norm_unref'][gO_k]) * n_k[gO_k] / voxel_volume

      if gO_sw[gO_k] > 0.0:
        dTStrans_norm [gO_k] = -GAS_CONST_KCAL * reference_density * 300.0 * gO_sw[gO_k] * np.log(gO_sw[gO_k]) #/ (gO_sw[gO_k] * reference_density) #* n_k[gO_k] / voxel_volume
      else:
        dTStrans_norm[gO_k] = 0.0

      if unhappy_only:

        if -(-Eww_norm[gO_k]+Esw_norm[gO_k])+(dTSorient_norm[gO_k]+dTStrans_norm [gO_k]) > 0.0:

          dTSorient_norm[gO_k]  = 0.0
          dTStrans_norm  [gO_k] = 0.0
          Esw_norm      [gO_k]  = 0.0
          Eww_norm      [gO_k]  = 0.0

      if happy_only:

        if -(-Eww_norm[gO_k]+Esw_norm[gO_k])+(dTSorient_norm[gO_k]+dTStrans_norm [gO_k]) < 0.0:

          dTSorient_norm[gO_k]  = 0.0
          dTStrans_norm  [gO_k] = 0.0
          Esw_norm      [gO_k]  = 0.0
          Eww_norm      [gO_k]  = 0.0

  structure_MS   = get_MS(path2molecule, radius_add=radius_add)

##MS
  volume_MS      = integrate_finite_tets(bins, crd, structure_MS, radius_add,  {}, is_volume_field=True)
  g_local_MS     = integrate_finite_tets(bins, crd, structure_MS, radius_add, grid_data=gO_sw) #/ volume_MS
  density_MS     = g_local_MS * reference_density
  n_total_MS     = density_MS * volume_MS

  Eww_self_MS    = 0.0
  if not Eww_matrix == None:
    #Eww_self_MS    = integrate_finite_tets(bins, crd, structure_MS, radius_add, grid_data=Eww_matrix, is_self_grid=True) / (g_local_MS * reference_density * volume_MS)
    #Eww_self_MS    = rough_integration_Eww (bins, crd, structure_MS, Eww_matrix, voxel_volume) / (g_local_MS * reference_density * volume_MS)
    Eww_self_MS    = rough_integration_Eww_smart (gist_data, structure_MS, Eww_matrix, voxel_volume)/(reference_density * g_local_MS * volume_MS)

  Eww_MS         = integrate_finite_tets(bins, crd, structure_MS, radius_add,  grid_data=Eww_norm) - Eww_self_MS #/n_total_MS

  Esw_MS         = - integrate_finite_tets(bins, crd, structure_MS, radius_add, grid_data=Esw_norm)       #/ volume_MS * n_total_MS 
  dTStrans_MS    = - integrate_finite_tets(bins, crd, structure_MS, radius_add, grid_data=dTStrans_norm)  #/ volume_MS * n_total_MS
  dTSorient_MS   = - integrate_finite_tets(bins, crd, structure_MS, radius_add, grid_data=dTSorient_norm) #/ volume_MS * n_total_MS
  delta_H_MS     = Esw_MS       + Eww_MS
  Tdelta_S_MS    = dTSorient_MS + dTStrans_MS
  delta_G_MS     = delta_H_MS   - Tdelta_S_MS


### rough integration ###
#MS
#we need a bounding box...
#  o           = np.array(crd[0,0,0])
#  xdim        = np.array(crd[bins[0]-1,0,0]) - o
#  ydim        = np.array(crd[0,bins[1]-1,0]) - o
#  zdim        = np.array(crd[0,0,bins[2]-1]) - o
#
#  bounding_min = np.array([min(structure_MS[:,0]), min(structure_MS[:,1]), min(structure_MS[:,2])]) - o
#  bounding_max = np.array([max(structure_MS[:,0]), max(structure_MS[:,1]), max(structure_MS[:,2])]) - o
#
#  min_x     = 0       if int(bounding_min[0]/xdim[0]*bins[0])-4 < 0       else int(bounding_min[0]/xdim[0]*bins[0])-4
#  min_y     = 0       if int(bounding_min[1]/ydim[1]*bins[1])-4 < 0       else int(bounding_min[1]/ydim[1]*bins[1])-4
#  min_z     = 0       if int(bounding_min[2]/zdim[2]*bins[2])-4 < 0       else int(bounding_min[2]/zdim[2]*bins[2])-4
#  max_x     = bins[0] if int(bounding_max[0]/xdim[0]*bins[0])+4 > bins[0] else int(bounding_max[0]/xdim[0]*bins[0])+4
#  max_y     = bins[1] if int(bounding_max[1]/ydim[1]*bins[1])+4 > bins[1] else int(bounding_max[1]/ydim[1]*bins[1])+4
#  max_z     = bins[2] if int(bounding_max[2]/zdim[2]*bins[2])+4 > bins[2] else int(bounding_max[2]/zdim[2]*bins[2])+4
#
#  tess_hull = Delaunay(structure_MS)
#  inside_coordinates = {}
#  for x in range(min_x, max_x):
#    for y in range(min_y, max_y):
#      for z in range(min_z, max_z):
#        p = crd[x,y,z]
#        if tess_hull.find_simplex(p)>=0:
#          inside_coordinates[x,y,z] = p
#
#  volume_MS      = rough_integration(bins, inside_coordinates, structure_MS, {}, voxel_volume)
#  g_local_MS     = rough_integration(bins, inside_coordinates, structure_MS, gO_sw, voxel_volume) / volume_MS
#  density_MS     = g_local_MS * reference_density
#  n_total_MS     = density_MS * volume_MS
#
#  #Eww_self_MS    = integrate_finite_tets(bins, crd, structure_MS, radius_add, grid_data=Eww_matrix, is_self_grid=True) / (g_local_MS * reference_density * volume_MS)
#  Eww_self_MS    = rough_integration_Eww (bins, inside_coordinates, structure_MS, Eww_matrix, voxel_volume) / (g_local_MS * reference_density * volume_MS)
#
#  Eww_MS         = 0.0
#  if Eww_matrix != {}:
#    Eww_MS         = 2.0*rough_integration(bins, inside_coordinates, structure_MS, Eww_norm)/volume_MS - 2.0*Eww_self_MS/volume_MS - reference_energy
#
#  Esw_MS         = rough_integration(bins, inside_coordinates, structure_MS, Esw_norm, voxel_volume)       / volume_MS
#  dTStrans_MS    = rough_integration(bins, inside_coordinates, structure_MS, dTStrans_norm, voxel_volume)  / volume_MS
#  dTSorient_MS   = rough_integration(bins, inside_coordinates, structure_MS, dTSorient_norm, voxel_volume) / volume_MS
#  delta_H_MS     = Esw_MS       + Eww_MS
#  Tdelta_S_MS    = dTSorient_MS + dTStrans_MS
#  delta_G_MS     = delta_H_MS   - Tdelta_S_MS

  if return_values:
    returner = OrderedDict()
    returner.update({'vol':volume_MS, 'g':g_local_MS, 'dens':density_MS, 'n_tot':n_total_MS, 'Eww_self':Eww_self_MS, 'Eww':Eww_MS, 'Esw':Esw_MS, \
            'dTStr':dTStrans_MS, 'dTSor':dTSorient_MS, 'dH':delta_H_MS, 'TdS':Tdelta_S_MS, 'dG':delta_G_MS})
    return returner

  warning   = ''
  if Eww_matrix=={}:
    warning = "Eww in volume is zero because there was no Eww-matrix."

  summary ='''
##  S u m m a r y  ##

general information:
%s
Ref Eww     : %6.3f
Ref dens    : %6.3f
Rad add     : %6.3f

# Molecular surface boundary #
total volume: %6.3f
total No wat: %6.3f
mean density: %6.3f
gO_average  : %6.3f

dEww        : %6.3f
dEww_self   : %6.3f
dEsw        : %6.3f
dH          : %6.3f

TdStrans    : %6.3f
TdSorient   : %6.3f
TdS         : %6.3f

dG          : %6.3f



''' %(warning, reference_energy, reference_density, radius_add, volume_MS, n_total_MS, density_MS, g_local_MS, Eww_MS, Eww_self_MS, Esw_MS, delta_H_MS, dTStrans_MS, dTSorient_MS, Tdelta_S_MS, delta_G_MS)

  print ("Writing output files...")
  o = open('%s_gist_cavity_analysis.dat' %output, "w")
  o.write(summary)
  o.close()
  print ("Done")

def rough_integration_sphere (bins, coordinate_data, radius, sphere_origin, grid_data, dV):

#we need a bounding box...
  o           = np.array(coordinate_data[0,0,0])
  xdim        = np.array(coordinate_data[bins[0]-1,0,0]) - o
  ydim        = np.array(coordinate_data[0,bins[1]-1,0]) - o
  zdim        = np.array(coordinate_data[0,0,bins[2]-1]) - o

  bounding_min = np.array(sphere_origin)-radius - o
  bounding_max = np.array(sphere_origin)+radius - o

  min_x     = 0       if int(bounding_min[0]/xdim[0]*bins[0])-4 < 0       else int(bounding_min[0]/xdim[0]*bins[0])-4
  min_y     = 0       if int(bounding_min[1]/ydim[1]*bins[1])-4 < 0       else int(bounding_min[1]/ydim[1]*bins[1])-4
  min_z     = 0       if int(bounding_min[2]/zdim[2]*bins[2])-4 < 0       else int(bounding_min[2]/zdim[2]*bins[2])-4
  max_x     = bins[0] if int(bounding_max[0]/xdim[0]*bins[0])+4 > bins[0] else int(bounding_max[0]/xdim[0]*bins[0])+4
  max_y     = bins[1] if int(bounding_max[1]/ydim[1]*bins[1])+4 > bins[1] else int(bounding_max[1]/ydim[1]*bins[1])+4
  max_z     = bins[2] if int(bounding_max[2]/zdim[2]*bins[2])+4 > bins[2] else int(bounding_max[2]/zdim[2]*bins[2])+4

  data = []
  for x in range(min_x, max_x):
    for y in range(min_y, max_y):
      for z in range(min_z, max_z):
        if calc_dist(coordinate_data[x,y,z], sphere_origin) <= radius:
          if grid_data == {}:
            data.append(1.0)
          else:
            data.append(grid_data[x,y,z])
  return (sum(data)*dV, np.sqrt(np.var(data))*dV)

def rough_integration (bins, coordinate_data, hull_object, grid_data, dV):

  data = 0.0
  for frac_crd, real_crd in coordinate_data.items():
    if grid_data == {}:
      data += dV
    if frac_crd in grid_data.keys():
      data += grid_data[frac_crd] * dV

  return data


def rough_integration_Eww_smart(gist_data, hull_object, path2Eww, dV):

  crd            = gist_data['crd']
  bins           = np.array(gist_data['bins'])

#we need a bounding box...
  o           = np.array(crd[0,0,0])
  xdim        = np.array(crd[bins[0]-1,0,0]) - o
  ydim        = np.array(crd[0,bins[1]-1,0]) - o
  zdim        = np.array(crd[0,0,bins[2]-1]) - o

  bounding_min = np.array([min(hull_object[:,0]), min(hull_object[:,1]), min(hull_object[:,2])]) - o
  bounding_max = np.array([max(hull_object[:,0]), max(hull_object[:,1]), max(hull_object[:,2])]) - o

  min_x     = 0       if int(bounding_min[0]/xdim[0]*bins[0])-4 < 0       else int(bounding_min[0]/xdim[0]*bins[0])-4
  min_y     = 0       if int(bounding_min[1]/ydim[1]*bins[1])-4 < 0       else int(bounding_min[1]/ydim[1]*bins[1])-4
  min_z     = 0       if int(bounding_min[2]/zdim[2]*bins[2])-4 < 0       else int(bounding_min[2]/zdim[2]*bins[2])-4
  max_x     = bins[0] if int(bounding_max[0]/xdim[0]*bins[0])+4 > bins[0] else int(bounding_max[0]/xdim[0]*bins[0])+4
  max_y     = bins[1] if int(bounding_max[1]/ydim[1]*bins[1])+4 > bins[1] else int(bounding_max[1]/ydim[1]*bins[1])+4
  max_z     = bins[2] if int(bounding_max[2]/zdim[2]*bins[2])+4 > bins[2] else int(bounding_max[2]/zdim[2]*bins[2])+4

  tess_hull = Delaunay(hull_object)
  inside_coordinates = {}
  for x in range(min_x, max_x):
    for y in range(min_y, max_y):
      for z in range(min_z, max_z):
        p = crd[x,y,z]
        if tess_hull.find_simplex(p)>=0:
          inside_coordinates[x,y,z] = p

  return crop_Eww_matrix(path2Eww, gist_data, crop_grid=inside_coordinates) * dV * len(inside_coordinates)


def rough_integration_Eww (bins, crd, hull_object, Eww_data, dV):

#we need a bounding box...
  o           = np.array(crd[0,0,0])
  xdim        = np.array(crd[bins[0]-1,0,0]) - o
  ydim        = np.array(crd[0,bins[1]-1,0]) - o
  zdim        = np.array(crd[0,0,bins[2]-1]) - o

  bounding_min = np.array([min(hull_object[:,0]), min(hull_object[:,1]), min(hull_object[:,2])]) - o
  bounding_max = np.array([max(hull_object[:,0]), max(hull_object[:,1]), max(hull_object[:,2])]) - o

  min_x     = 0       if int(bounding_min[0]/xdim[0]*bins[0])-4 < 0       else int(bounding_min[0]/xdim[0]*bins[0])-4
  min_y     = 0       if int(bounding_min[1]/ydim[1]*bins[1])-4 < 0       else int(bounding_min[1]/ydim[1]*bins[1])-4
  min_z     = 0       if int(bounding_min[2]/zdim[2]*bins[2])-4 < 0       else int(bounding_min[2]/zdim[2]*bins[2])-4
  max_x     = bins[0] if int(bounding_max[0]/xdim[0]*bins[0])+4 > bins[0] else int(bounding_max[0]/xdim[0]*bins[0])+4
  max_y     = bins[1] if int(bounding_max[1]/ydim[1]*bins[1])+4 > bins[1] else int(bounding_max[1]/ydim[1]*bins[1])+4
  max_z     = bins[2] if int(bounding_max[2]/zdim[2]*bins[2])+4 > bins[2] else int(bounding_max[2]/zdim[2]*bins[2])+4

  tess_hull = Delaunay(hull_object)
  inside_coordinates = {}
  for x in range(min_x, max_x):
    for y in range(min_y, max_y):
      for z in range(min_z, max_z):
        p = crd[x,y,z]
        if tess_hull.find_simplex(p)>=0:
          inside_coordinates[x,y,z] = p

  Eww = 0.0
  for frac_crd1, real_crd1 in inside_coordinates.items():
    for frac_crd2, real_crd2 in inside_coordinates.items():
      if frac_crd2 in Eww_data[frac_crd1].keys():
        Eww += Eww_data[frac_crd1][frac_crd2] * dV
  return Eww*dV


def integrate_finite_tets(bins, coordinate_data, hull_object,  radius_add, grid_data={}, is_self_grid=False, is_volume_field=False):
#sampling_rate gibt an ueber wieviele punkte auf einer tetraederachse(baryzentric!) integriert wird, alpha ist der filter-faktor fuer alpha-shapes
#  self_volume = 0.0
#  if is_self_grid:
#    self_volume = integrate_finite_tets(bins, coordinate_data, hull_object, grid_data={}, is_self_grid=False,is_volume_field=True)

  sampling_size = 0.5
  #alpha         = 5.0+radius_add**2

#we need a bounding box...
  o           = np.array(coordinate_data[0,0,0])
  xdim        = np.array(coordinate_data[bins[0]-1,0,0]) - o
  ydim        = np.array(coordinate_data[0,bins[1]-1,0]) - o
  zdim        = np.array(coordinate_data[0,0,bins[2]-1]) - o

  bounding_min = np.array([min(hull_object[:,0]), min(hull_object[:,1]), min(hull_object[:,2])]) - o
  bounding_max = np.array([max(hull_object[:,0]), max(hull_object[:,1]), max(hull_object[:,2])]) - o

  min_x     = 0       if int(bounding_min[0]/xdim[0]*bins[0])-4 < 0       else int(bounding_min[0]/xdim[0]*bins[0])-4
  min_y     = 0       if int(bounding_min[1]/ydim[1]*bins[1])-4 < 0       else int(bounding_min[1]/ydim[1]*bins[1])-4
  min_z     = 0       if int(bounding_min[2]/zdim[2]*bins[2])-4 < 0       else int(bounding_min[2]/zdim[2]*bins[2])-4
  max_x     = bins[0] if int(bounding_max[0]/xdim[0]*bins[0])+4 > bins[0] else int(bounding_max[0]/xdim[0]*bins[0])+4
  max_y     = bins[1] if int(bounding_max[1]/ydim[1]*bins[1])+4 > bins[1] else int(bounding_max[1]/ydim[1]*bins[1])+4
  max_z     = bins[2] if int(bounding_max[2]/zdim[2]*bins[2])+4 > bins[2] else int(bounding_max[2]/zdim[2]*bins[2])+4

#save some memory...
  del o, xdim, ydim, zdim, bounding_min, bounding_max

#Then we build the sampling grid...
  xx = np.array(map(lambda x: coordinate_data[x,0,0][0], range(min_x, max_x)))
  yy = np.array(map(lambda y: coordinate_data[0,y,0][1], range(min_y, max_y)))
  zz = np.array(map(lambda z: coordinate_data[0,0,z][2], range(min_z, max_z)))

  def grid_the_values(x,y,z):
    if is_volume_field:
      return 1.0
    if is_self_grid and (x,y,z) in grid_data:
      return integrate_finite_tets(bins, coordinate_data, hull_object, radius_add, grid_data=grid_data[x,y,z], is_self_grid=False)
    if not (x,y,z) in grid_data:
      return 0.0    
    return grid_data[x,y,z]

  grid_data_vectorized             = np.vectorize(grid_the_values)  
  data                             = grid_data_vectorized(*np.meshgrid(np.array(range(min_x,max_x),dtype=float), np.array(range(min_y,max_y), dtype=float), \
                                                                                     np.array(range(min_z,max_z),dtype=float), indexing='ij', sparse=True))
  value_field_interpolation_linear = scipy.interpolate.RegularGridInterpolator((xx,yy,zz), data, bounds_error=False, fill_value=None)

  tess_hull   = Delaunay(hull_object)
  tets        = tess_hull.points[tess_hull.simplices]
  tets_volume = np.abs(np.einsum('ij,ij->i', tets[:,0]-tets[:,3], np.cross(tets[:,1]-tets[:,3], tets[:,2]-tets[:,3]))) / 6
#calculate alpha shapes. Alpha is unit angstrom
  a = np.linalg.norm(tets[:,0] - tets[:,1], axis=1)
  b = np.linalg.norm(tets[:,0] - tets[:,2], axis=1)
  c = np.linalg.norm(tets[:,0] - tets[:,3], axis=1)
  p = np.linalg.norm(tets[:,3] - tets[:,2], axis=1)
  q = np.linalg.norm(tets[:,3] - tets[:,1], axis=1)
  r = np.linalg.norm(tets[:,2] - tets[:,1], axis=1)

  radii = np.sqrt(abs( (a*p+b*q+c*r) * (-a*p+b*q+c*r) * (a*p-b*q+c*r) * (a*p+b*q-c*r)) ) / (24.0*tets_volume)

  tet_is_in_alpha = np.ones(len(tets),dtype=bool)
#  for tet_index in range(len(tets)):
#
#    if radii[tet_index] > alpha:
#      tet_is_in_alpha[tet_index] = False

#  if is_volume_field:
#    print "Excluded %d from %d simplices with alpha=%4.2f." %(len(np.where(tet_is_in_alpha==False)[0]), len(tet_is_in_alpha), alpha)
#    print "Now using sampling_size=%4.2f for integration." %sampling_size

#  temp = {}
  Integral = 0.0
  for index in range(len(tets)):

    if not tet_is_in_alpha[index]:
      continue

    if np.isnan(tess_hull.transform[index]).all():
      if is_volume_field: print "Warning: weird simplex detected."
      continue

    b0_sampling_rate = int(distance.euclidean(tess_hull.points[tess_hull.simplices[index][0]], tess_hull.points[tess_hull.simplices[index][3]])/sampling_size)+2
    b1_sampling_rate = int(distance.euclidean(tess_hull.points[tess_hull.simplices[index][1]], tess_hull.points[tess_hull.simplices[index][3]])/sampling_size)+2
    b2_sampling_rate = int(distance.euclidean(tess_hull.points[tess_hull.simplices[index][2]], tess_hull.points[tess_hull.simplices[index][3]])/sampling_size)+2

    b0_sampling           = np.linspace(0.0,1.0,b0_sampling_rate, dtype=float)
    b1_sampling           = np.linspace(0.0,1.0,b1_sampling_rate, dtype=float)
    b2_sampling           = np.linspace(0.0,1.0,b2_sampling_rate, dtype=float)
    baryzentric_grid_real = np.meshgrid(b0_sampling, b1_sampling, b2_sampling, indexing='ij')

    xyz          = np.einsum('ij,jhkl->ihkl', np.linalg.inv(tess_hull.transform[index,:3]), baryzentric_grid_real)
    x            = xyz[0] + tess_hull.transform[index,3][0]
    y            = xyz[1] + tess_hull.transform[index,3][1]
    z            = xyz[2] + tess_hull.transform[index,3][2]
    sampled_data = value_field_interpolation_linear((x,y,z))

#   pdb-writing for debugging...
#    for b0 in range(b0_sampling_rate):
#      b1_limit = b1_sampling_rate - b0
#      for b1 in range(b1_limit):
#        b2_limit = b2_sampling_rate - b0 - b1
#        for b2 in range(b2_limit):
#          temp[x[b0,b1,b2], y[b0,b1,b2], z[b0,b1,b2]] = 1.0


    Integral_b0  = np.zeros(b0_sampling_rate)
    Integral_b1  = np.zeros(b1_sampling_rate)
    for b0 in range(b0_sampling_rate):
      b1_limit = b1_sampling_rate - b0
      for b1 in range(b1_limit):
        b2_limit = b2_sampling_rate - b0 - b1
        Integral_b1[b1] = np.trapz(sampled_data[b0,b1,0:b2_limit], b2_sampling[0:b2_limit])
      Integral_b0[b0]=np.trapz(Integral_b1[0:b1_limit], b1_sampling[0:b1_limit])
    Integral += np.trapz(Integral_b0, b0_sampling)

#  write_pdb(temp,"tetscrds_factor_5.0")
#  exit()
  return Integral

def gist_grid_analysis (gist_data, Eww_matrix, calculation_density, reference_density, reference_energy):
#this method returns average values over a whole GIST-box. Useful for obtaining Bulk-references

  #dTStrans_norm  = gist_data['dTStrans_norm']
  dTSorient_norm = gist_data['dTSorient_norm']
  Esw_norm       = gist_data['Esw_norm']
  Eww_norm       = gist_data['Eww_norm_unref']

  gO_GIST        = gist_data['gO']
  gO_sw          = {}
#  gO_sw_k_crds   = np.array(gO_sw_k.values())[:,0]
#  gO_sw_k_value  = np.array(gO_sw_k.values())[:,1]

  bins           = gist_data['bins']
  n              = gist_data['n']
  dim            = gist_data['dim']

#dV is the voxel_volume
  dV      = (dim[0] / n[0]) * (dim[1] / n[1]) * (dim[2] / n[2])

  dTStrans  = []
  dTSorient = []
  Eww       = []
  Esw       = []
  density   = []
  dH        = []
  TdS       = []
  dG        = []
  norm_fac  = []

#Integrate over all voxels k
  for x in range(0, bins[0]):
    for y in range(0, bins[1]):
      for z in range(0, bins[2]):

#integrate interaction energies with all other voxels l
#        Eww_self_crds   = np.array(Eww_matrix[x,y,z].items())[,0]
#        Eww_self_values = np.array(Eww_matrix[x,y,z].items())[,1]
#        gO_sw_l         = gO_sw_k_crds[np.in1d(gO_sw_k_crds, Eww_self_crds)]

#correct gO to reference density
        density_k        = gO_GIST[x,y,z]               * calculation_density
        gO_sw            = density_k / reference_density
        Eww_k_self       = 0.0
        if (x,y,z) in Eww_matrix:
          # "Integrate"
          Eww_k_self   = 2.0 * sum(Eww_matrix[x,y,z].values())
          # "transform to per water normalization"
          Eww_k_self   /= gO_sw * reference_density

        Eww_k            = 2.0*Eww_norm[x,y,z]   * dV
        Esw_k            = Esw_norm[x,y,z]       * dV
        dTSorient_k      = dTSorient_norm[x,y,z] * dV

        dTStrans_dens_k  = - GAS_CONST_KCAL * reference_density * 300.0 * gO_sw * np.log(gO_sw)
        dTStrans_k       = dTStrans_dens_k / (gO_sw * reference_density) * dV

        Eww.append(Eww_k - Eww_k_self - reference_energy)
        Esw.append(Esw_k)
        dTStrans.append(dTStrans_k)
        dTSorient.append(dTSorient_k)
        density.append(density_k*dV)
        dH.append(Eww_k+Esw_k)
        TdS.append(dTStrans_k+dTSorient_k)
        dG.append(Eww_k+Esw_k-(dTStrans_k+dTSorient_k))
        norm_fac.append(dV)

  voxel_counts = float(bins[0] * bins[1] * bins[2])
  grid_volume  = dV * voxel_counts
  Eww_vol      = (sum(Eww)      /sum(norm_fac), np.sqrt(np.var(Eww))      /sum(norm_fac)*10000)
  Esw_vol      = (sum(Esw)      /sum(norm_fac), np.sqrt(np.var(Esw))      /sum(norm_fac)*10000)
  TdS_tr_vol   = (sum(dTStrans) /sum(norm_fac), np.sqrt(np.var(dTStrans)) /sum(norm_fac)*10000)
  TdS_or_vol   = (sum(dTSorient)/sum(norm_fac), np.sqrt(np.var(dTSorient))/sum(norm_fac)*10000)
  dens_vol     = (sum(density)  /sum(norm_fac), np.sqrt(np.var(density))  /sum(norm_fac)*10000)
  dH_vol       = (sum(dH)       /sum(norm_fac), np.sqrt(np.var(dH))       /sum(norm_fac)*10000) 
  TdS_vol      = (sum(TdS)      /sum(norm_fac), np.sqrt(np.var(TdS))      /sum(norm_fac)*10000)
  dG_vol       = (sum(dG)       /sum(norm_fac), np.sqrt(np.var(dG))       /sum(norm_fac)*10000)
  n_wat_vol    = (sum(density)                , np.sqrt(np.var(density))                *10000)

  print '\n'
  print '#### Normalized properties on grid ###'
  print '                   Eww[kcal/mol] Esw[kcal/mol] TdS_tr[kcal/mol] TdS_or[kcal/mol] dens[wat-mol/A^3] n_wat       dH[kcal/mol]  TdS[kcal/mol] dG[kcal/mol]'
  print 'Value            :  %6.5f      %6.5f     %6.5f          %6.5f          %6.5f           %6.5f   %6.5f      %6.5f       %6.5f       ' %(Eww_vol[0], \
   Esw_vol[0], TdS_tr_vol[0],  TdS_or_vol[0], dens_vol[0], n_wat_vol[0], dH_vol[0], TdS_vol[0], dG_vol[0])
  print 'Error [10^-4] +/-:    %6.5f      %6.5f     %6.5f          %6.5f          %6.5f           %6.5f     %6.5f       %6.5f       %6.5f       ' %(Eww_vol[1], \
   Esw_vol[1], TdS_tr_vol[1], TdS_or_vol[1], dens_vol[1], n_wat_vol[1], dH_vol[1], TdS_vol[1], dG_vol[1])

  print 'Total grid volume:   %8.2f A^3'  %grid_volume


def population_analysis_michael (water_positions, grid_data, radius, output):

  bins            = np.array(grid_data['bins'])
  n               = np.array(grid_data['n'])
  dim             = np.array(grid_data['dim'])
  origin          = np.array(grid_data['origin'])
  alpha           = np.array(grid_data['alpha'])
  beta            = np.array(grid_data['beta'])
  gamma           = np.array(grid_data['gamma'])

  population_data = grid_data['map']
  propensity_data = {}
  coordinates     = OrderedDict()
  voxel_volume   = (dim[0] / n[0]) * (dim[1] / n[1]) * (dim[2] / n[2])
  
  cos_a = ( math.cos(math.pi*beta/180) * math.cos(math.pi*gamma/180) - math.cos(math.pi*alpha/180) ) / ( math.sin(math.pi*beta/180) *  math.sin(math.pi*gamma/180))
  sin_a = math.sqrt(1.0 - cos_a**2)

  for key, value in population_data.items():

    X_coord = origin[0] + dim[0] / n[0] * key[0] + dim[1] / n[1] * math.cos(math.pi*gamma/180) * key[1] + dim[2] / n[2] * math.cos(math.pi*beta/180) * key[2]
    Y_coord = origin[1] + 0                      + dim[1] / n[1] * math.sin(math.pi*gamma/180) * key[1] - dim[2] / n[2] * math.sin(math.pi*beta/180) * cos_a * key[2]
    Z_coord = origin[2] + 0                      + 0                                                    + dim[2] / n[2] * math.sin(math.pi*beta/180) * sin_a * key[2]
    coordinates[key[0], key[1], key[2]] = [X_coord, Y_coord, Z_coord]

  for crd, value in population_data.items():
    propensity_data[crd[0], crd[1], crd[2]] = value/voxel_volume

  output_pop_data  = OrderedDict()
  output_prop_data = OrderedDict()
  dict_indexing    = OrderedDict()

  x_bounding_limit = radius / dim[0]
  y_bounding_limit = radius / dim[1]
  z_bounding_limit = radius / dim[2]

  for index, water_position in enumerate(water_positions.keys()):

    box_checker = is_in_box(water_position, grid_data, return_coeffs=True)[1]

    if not (x_bounding_limit < box_checker[0] < 1.0-x_bounding_limit and y_bounding_limit < box_checker[1] < 1.0-y_bounding_limit and z_bounding_limit < box_checker[2] < 1.0-z_bounding_limit):
      print "Warning: Hydration site %d out of boundaries. You'd better lower integration radius, use larger grid or just deal with it." %(index+1)
      continue

    volume                           = integrate_finite_spheres(bins, coordinates, water_position, radius, grid_data={}, is_volume_field=True)[0]
    output_pop_data[water_position]  = integrate_finite_spheres(bins, coordinates, water_position, radius, grid_data=population_data)[0] / volume
    output_prop_data[water_position] = integrate_finite_spheres(bins, coordinates, water_position, radius, grid_data=propensity_data)[0]
    dict_indexing[water_position]    = index+1

  write_pdb(output_pop_data ,"pop_" +output,dict_indexing)
  #write_pdb(output_prop_data,"prop_"+output,dict_indexing)


def gist_water_analysis (water_positions, gist_data, radius, Eww_matrix, reference_energy, calculation_density, reference_density, output, return_it=False, radius_start=0.0, gist_data_apo=None):

#  if gist_data_apo != None:
#    apo_data      = gist_water_analysis(water_positions, gist_data_apo, radius, Eww_matrix, reference_energy, calculation_density, reference_density, output, return_it=True, radius_start=0.0)
#    n_tot_apo     = apo_data['n_tot']
#    dens_apo      = apo_data['dens']
#    dH_apo        = apo_data['dH']
#    TdS_apo       = apo_data['TdS']
#    dG_apo        = apo_data['dG']
#    dEww_apo      = apo_data['dEww']
#    dEww_self_apo = apo_data['dEww_self']
#    dEsw_apo      = apo_data['dEsw']
#    dTStra_apo    = apo_data['dTStra']
#    dTSor_apo     = apo_data['dTSor']

  bins           = np.array(gist_data['bins'])
  n              = np.array(gist_data['n'])
  dim            = np.array(gist_data['dim'])

  gO_GIST        = gist_data['gO']

  coordinates    = gist_data['crd']

  voxel_volume   = (dim[0] / n[0]) * (dim[1] / n[1]) * (dim[2] / n[2])

#correct gO, calculate Entropy_trans and correct with reference energy
  gO_sw          = {}
  dTStrans_norm  = {}
  dTSorient_norm = {}
  Esw_norm       = {}
  Eww_norm       = {}
  n_k            = {}

  for gO_k, gO_k_value in gO_GIST.items():
    gO_sw[gO_k]          = gO_k_value  * calculation_density / reference_density
    n_k  [gO_k]          = gO_sw[gO_k] * reference_density * voxel_volume
    dTSorient_norm[gO_k] = gist_data['dTSorient_norm'][gO_k] * n_k[gO_k] / voxel_volume
    Esw_norm      [gO_k] = gist_data['Esw_norm'][gO_k] * n_k[gO_k] / voxel_volume
    Eww_norm      [gO_k] = (reference_energy - 2.0*gist_data['Eww_norm_unref'][gO_k]) * n_k[gO_k] / voxel_volume

    if gO_sw[gO_k] > 0.0:
      dTStrans_norm [gO_k] = -GAS_CONST_KCAL * reference_density * 300.0 * gO_sw[gO_k] * np.log(gO_sw[gO_k]) #/ (gO_sw[gO_k] * reference_density) #* n_k[gO_k] / voxel_volume
    else:
      dTStrans_norm[gO_k] = 0.0

  water_position_n_average  = OrderedDict()
  water_position_dens       = OrderedDict()
  water_positions_Eww       = OrderedDict()
  water_positions_Eww_self  = OrderedDict()
  water_positions_Esw       = OrderedDict()
  water_positions_dTStrans  = OrderedDict()
  water_positions_dTSorient = OrderedDict()
  water_positions_delta_H   = OrderedDict()
  water_positions_delta_TS  = OrderedDict()
  water_positions_delta_G   = OrderedDict()
  dict_indexing             = OrderedDict()
  header = 'Site#   vol[A^3] dEww_self   dEww   dEsw   TdStrans  dSorient    dH       TdS      dG      n_avr    dens[wat/A^3]  g_O    FWHM[A]\n'
  data   = ''
  data_error = ''

  x_bounding_limit = (radius+3.0*dim[0]/n[0]) / dim[0]
  y_bounding_limit = (radius+3.0*dim[1]/n[1]) / dim[1]
  z_bounding_limit = (radius+3.0*dim[2]/n[2]) / dim[2]

  print "Calculating water energies..."
  for index, water_position in enumerate(water_positions.keys()):

    box_checker = is_in_box(water_position, gist_data, return_coeffs=True)[1]

    if not (x_bounding_limit < box_checker[0] < 1.0-x_bounding_limit and y_bounding_limit < box_checker[1] < 1.0-y_bounding_limit and z_bounding_limit < box_checker[2] < 1.0-z_bounding_limit):
      data       +=  ' %3d   Not on grid or too near to boundary.\n' %(index+1)
      data_error +=  ' %3d   Not on grid or too near to boundary.\n' %(index+1)
      continue
    #First g_local calculation is just a dummy for the calculation of FWHM
    g_local_temp, g_local_half_maximum,  g_local_FWHM  = integrate_finite_spheres(bins, coordinates, water_position, radius+3.0*dim[0]/n[0], radius_start=radius_start, grid_data=gO_sw, FWHM=True)

    volume    = integrate_finite_spheres(bins, coordinates, water_position, radius, radius_start=radius_start, grid_data={}, is_volume_field=True)[0]#rough_integration_sphere (bins, coordinates, radius, water_position, {}, voxel_volume)
    g_local   = integrate_finite_spheres(bins, coordinates, water_position, radius, radius_start=radius_start, grid_data=gO_sw)[0]/volume#rough_integration_sphere (bins, coordinates, radius, water_position, gO_sw, voxel_volume)[0]/volume
    n_average = integrate_finite_spheres(bins, coordinates, water_position, radius, radius_start=radius_start, grid_data=n_k)[0]/volume#rough_integration_sphere (bins, coordinates, radius, water_position, n_k, voxel_volume)[0]/volume
    density   = n_average / volume

    Eww_self       = 0.0
    Eww_self_error = 0.0

    if not Eww_matrix=={}:
      Eww_self        = integrate_finite_spheres(bins, coordinates, water_position, radius, radius_start=radius_start,  grid_data=Eww_matrix, is_self_grid=True)[0] * n_average
      Eww_self_error  = integrate_finite_spheres(bins, coordinates, water_position, radius+3.0*dim[0]/n[0], radius_start=radius_start,  grid_data=Eww_matrix, is_self_grid=True, FWHM=True)[1] * n_average * volume

    #Eww        = (rough_integration_sphere (bins, coordinates, radius, water_position, Eww_norm, voxel_volume)[0] - Eww_self) #* n_average
    #Esw        = -rough_integration_sphere (bins, coordinates, radius, water_position, Esw_norm, voxel_volume)[0] #* n_average
    #dTStrans   = -rough_integration_sphere (bins, coordinates, radius, water_position, dTStrans_norm, voxel_volume)[0]  #* n_average
    #dTSorient  = -rough_integration_sphere (bins, coordinates, radius, water_position, dTSorient_norm, voxel_volume)[0] #* n_average

    Eww        = (integrate_finite_spheres(bins, coordinates, water_position, radius, radius_start=radius_start, grid_data=Eww_norm)[0]  - Eww_self)
    Esw        = - integrate_finite_spheres(bins, coordinates, water_position, radius, radius_start=radius_start, grid_data=Esw_norm)[0]
    dTStrans   = - integrate_finite_spheres(bins, coordinates, water_position, radius, radius_start=radius_start, grid_data=dTStrans_norm)[0]
    dTSorient  = - integrate_finite_spheres(bins, coordinates, water_position, radius, radius_start=radius_start, grid_data=dTSorient_norm)[0]

    Eww_error       = integrate_finite_spheres(bins, coordinates, water_position, radius+3.0*dim[0]/n[0], radius_start=radius_start, grid_data=Eww_norm, FWHM=True)[1]*volume
    Esw_error       = integrate_finite_spheres(bins, coordinates, water_position, radius+3.0*dim[0]/n[0], radius_start=radius_start, grid_data=Esw_norm, FWHM=True)[1]*volume
    dTStrans_error  = integrate_finite_spheres(bins, coordinates, water_position, radius+3.0*dim[0]/n[0], radius_start=radius_start, grid_data=dTStrans_norm, FWHM=True)[1]*volume
    dTSorient_error = integrate_finite_spheres(bins, coordinates, water_position, radius+3.0*dim[0]/n[0], radius_start=radius_start, grid_data=dTSorient_norm, FWHM=True)[1]*volume

    delta_H    = Esw + Eww
    Tdelta_S   = dTSorient + dTStrans
    delta_G    = delta_H - Tdelta_S
    water_positions_delta_H  [water_position]  = delta_H
    water_positions_delta_TS [water_position]  = Tdelta_S
    water_positions_delta_G  [water_position]  = delta_G
    water_positions_Eww      [water_position]  = Eww
    water_positions_Eww_self [water_position]  = Eww_self
    water_positions_Esw      [water_position]  = Esw
    water_positions_dTStrans [water_position]  = dTStrans
    water_positions_dTSorient[water_position]  = dTSorient
    water_position_n_average [water_position]  = n_average
    water_position_dens      [water_position]  = density
    dict_indexing            [water_position]  = index+1

    data += ' %3d    %6.3f    %6.3f   %6.3f  %6.3f   %6.3f    %6.3f   %6.3f   %6.3f   %6.3f   %6.3f    %6.3f    %6.3f    %6.3f\n'\
    %(index+1, volume, Eww_self, Eww, Esw, dTStrans, dTSorient, delta_H, Tdelta_S, delta_G, n_average, density, g_local, g_local_FWHM)

    data_error += ' %3d    %6.3s    %6.3f   %6.3f  %6.3f   %6.3f    %6.3f   %6.3s   %6.3s   %6.3s   %6.3s    %6.3s    %6.3f    %6.3f\n'\
    %(index+1, "NA", Eww_self_error, Eww_error, Esw_error, dTStrans_error, dTSorient_error, "NA", "NA", "NA", "NA", "NA", g_local_half_maximum*volume, g_local_FWHM)

  if return_it:
    returner = OrderedDict()
    returner.update({'n_avg':water_position_n_average, \
                     'dens':water_position_dens, \
                     'dH':water_positions_delta_H, \
                     'TdS':water_positions_delta_TS, \
                     'dG':water_positions_delta_G, \
                     'dEww':water_positions_Eww, \
                     'dEww_self':water_positions_Eww_self, \
                     'dEsw':water_positions_Esw, \
                     'dTStra':water_positions_dTStrans, \
                     'dTSor':water_positions_dTSorient, \
                     })

    return returner

  print ("Writing output files...")
  o = open('%s_gist_water_analysis.dat' %output, "w")
  o.write(header+data)
  o.close()
  o = open('%s_gist_water_analysis_errorestimation.dat' %output, "w")
  o.write(header+data_error)
  o.close()
  write_pdb(water_positions_delta_H  ,'%s_gist_water_analysis_dH'    %output, dict_indexing)
  write_pdb(water_positions_delta_TS ,'%s_gist_water_analysis_TdS'   %output, dict_indexing)
  write_pdb(water_positions_delta_G  ,'%s_gist_water_analysis_dG'    %output, dict_indexing)
  write_pdb(water_positions_Eww      ,'%s_gist_water_analysis_Eww'   %output, dict_indexing)
  write_pdb(water_positions_Esw      ,'%s_gist_water_analysis_Esw'   %output, dict_indexing)
  write_pdb(water_positions_dTStrans ,'%s_gist_water_analysis_dTStr' %output, dict_indexing)
  write_pdb(water_positions_dTSorient,'%s_gist_water_analysis_dTSor' %output, dict_indexing)
  print ("Done")

def integrate_finite_spheres(bins, coordinate_data, sphere_origin, radius, radius_start=0.0, grid_data={}, is_self_grid=False, is_volume_field=False, FWHM=False):
#if is_self_grid, grid_data must be dict[x,y,z] = {dict[x,y,z]:float(value)}
#we need a bounding box...
#  self_volume = 0.0
#  if is_self_grid:
#    self_volume = integrate_finite_spheres(bins, coordinate_data, sphere_origin, radius, radius_start, {}, is_self_grid=False,is_volume_field=True)

  sampling_grid_size = 0.1

  o            = np.array(coordinate_data[0,0,0])
  xdim         = np.array(coordinate_data[bins[0]-1,0,0]) - o
  ydim         = np.array(coordinate_data[0,bins[1]-1,0]) - o
  zdim         = np.array(coordinate_data[0,0,bins[2]-1]) - o

  bounding_min = np.array(sphere_origin)-radius - o
  bounding_max = np.array(sphere_origin)+radius - o

  min_x     = 0       if int(bounding_min[0]/xdim[0]*bins[0])-10 < 0       else int(bounding_min[0]/xdim[0]*bins[0])-10
  min_y     = 0       if int(bounding_min[1]/ydim[1]*bins[1])-10 < 0       else int(bounding_min[1]/ydim[1]*bins[1])-10
  min_z     = 0       if int(bounding_min[2]/zdim[2]*bins[2])-10 < 0       else int(bounding_min[2]/zdim[2]*bins[2])-10
  max_x     = bins[0] if int(bounding_max[0]/xdim[0]*bins[0])+10 > bins[0] else int(bounding_max[0]/xdim[0]*bins[0])+10
  max_y     = bins[1] if int(bounding_max[1]/ydim[1]*bins[1])+10 > bins[1] else int(bounding_max[1]/ydim[1]*bins[1])+10
  max_z     = bins[2] if int(bounding_max[2]/zdim[2]*bins[2])+10 > bins[2] else int(bounding_max[2]/zdim[2]*bins[2])+10

  sampling_rate = int(radius/sampling_grid_size)+1

#Then we build the sampling grid...
  xx = np.array(map(lambda x: coordinate_data[x,0,0][0], range(min_x, max_x)))
  yy = np.array(map(lambda y: coordinate_data[0,y,0][1], range(min_y, max_y)))
  zz = np.array(map(lambda z: coordinate_data[0,0,z][2], range(min_z, max_z)))

  def grid_the_values(x,y,z):
    if is_volume_field:
      return 1.0
    if is_self_grid and (x,y,z) in grid_data:
      return integrate_finite_spheres(bins, coordinate_data, sphere_origin, radius, grid_data=grid_data[x,y,z], is_self_grid=False)
    if not (x,y,z) in grid_data:
      return 0.0
    return grid_data[x,y,z]

  grid_data_vectorized             = np.vectorize(grid_the_values)  
  data                             = grid_data_vectorized(*np.meshgrid(np.array(range(min_x,max_x),dtype=float), np.array(range(min_y,max_y), dtype=float), \
                                                                                     np.array(range(min_z,max_z),dtype=float), indexing='ij', sparse=True))
  value_field_interpolation_linear = scipy.interpolate.RegularGridInterpolator((xx,yy,zz), data, bounds_error=False, fill_value=None)

  def sphere_to_xyz (r, theta, phi):
    x, y, z = r * np.sin(theta) * np.cos(phi) + sphere_origin[0], r * np.sin(theta) * np.sin(phi) + sphere_origin[1],  r * np.cos(theta) + sphere_origin[2]
    return value_field_interpolation_linear((x,y,z)) * r**2 * np.sin(theta)

  phi_sampling   = np.linspace(0.0, 2*np.pi, sampling_rate)
  theta_sampling = np.linspace(0.0, np.pi,   sampling_rate)
  r_sampling     = np.linspace(radius_start, radius,  sampling_rate)

  sampled_data = sphere_to_xyz(*np.meshgrid(r_sampling, theta_sampling, phi_sampling, indexing='ij', sparse=True))
  Integral_r     = np.zeros(sampling_rate)
  Integral_theta = np.zeros(sampling_rate)

  for r in range(sampling_rate):

    for theta in range(sampling_rate):

      Integral_theta[theta]=np.trapz(sampled_data[r,theta,:], phi_sampling)

    Integral_r[r]=np.trapz(Integral_theta, theta_sampling)

  result = np.trapz(Integral_r, r_sampling)

#calculate FWHM
  if not FWHM:
    return (result, None, None)

  big_radius    = 0.5*np.sqrt(3.0)*(max_x-min_x)*(xdim[0]/bins[0])
  sampling_rate = int(big_radius/sampling_grid_size)+1

  phi_sampling   = np.linspace(0.0, 2*np.pi, sampling_rate)
  theta_sampling = np.linspace(0.0, np.pi,   sampling_rate)
  r_sampling     = np.linspace(0.0, big_radius,  sampling_rate)

  sampled_data   = sphere_to_xyz(*np.meshgrid(r_sampling, theta_sampling, phi_sampling, indexing='ij', sparse=True))
  Integral_r     = np.zeros(sampling_rate)
  Integral_theta = np.zeros(sampling_rate)

  for r in range(sampling_rate):

    for theta in range(sampling_rate):

      Integral_theta[theta]=np.trapz(sampled_data[r,theta,:], phi_sampling)

    Integral_r[r]=np.trapz(Integral_theta, theta_sampling)

    if r>3:

      Integral_last=np.trapz(Integral_r[r-2:r], r_sampling[r-2:r])     / (4.0/3.0 * np.pi * (((r-1)*sampling_grid_size)**3-((r-2)*sampling_grid_size)**3))
      Integral_new =np.trapz(Integral_r[r-1:r+1], r_sampling[r-1:r+1]) / (4.0/3.0 * np.pi * ((r*sampling_grid_size)**3-((r-1)*sampling_grid_size)**3))

      if -0.1 < (Integral_new-Integral_last)/sampling_grid_size < 0.1:

        maximum      = np.trapz(Integral_r[1:3], r_sampling[1:3]) / (4.0/3.0 * np.pi * ((2*sampling_grid_size)**3-(sampling_grid_size)**3))
        half_maximum = maximum - 0.5*(maximum-Integral_new)

        for r1 in range(1, sampling_rate):

          value = np.trapz(Integral_r[r1:r1+2], r_sampling[r1:r+2])     / (4.0/3.0 * np.pi * (((r1+1)*sampling_grid_size)**3-(r1*sampling_grid_size)**3))

          if value<half_maximum:

            return (result, half_maximum, (r1+0.5)*sampling_grid_size)

    

#via gauss quadrature, not stable...
  #return scipy.integrate.tplquad(sphere_to_xyz, 0, radius, lambda x: 0.0, lambda x: np.pi, lambda x,y: 0.0, lambda x,y: 2*np.pi)[0]


def density_vs_time (traj, prmtop, h_site, site_radius, name, np, start_init, stop_init, skip):

  site_occupancy = OrderedDict()

  if np > 1:

    O_coord = MDAnalysis.Universe(h_site).coord
    h_site_count = len(O_coord)
    MDtraj = MDAnalysis.Universe(prmtop, traj)
    stop_init = MDtraj.trajectory[stop_init].frame
    del O_coord, MDtraj #saving some memory...
    length = stop_init - start_init
    processes = []
  
    for i in range(np):

      start = start_init + length / np * i
      stop = start_init + length / np * (i + 1)
      p = subprocess.Popen([sys.executable, path2script, '-j', 'dens_v_t', '-pa', '%s' %prmtop, '-tj', '%s' %traj, '-hi', '%s' %h_site, '-hr', '%s' %site_radius, '-np', '1', '-o', '%s' %name, '-start', '%s' %str(start), '-stop', '%s' %str(stop), '-skip', '%s' %str(skip)], shell = False)
      processes.append(p)

    for i in range(np):

      processes[i].wait()
      f  = open(".tmp_%s.dat" %str(start_init + length / np * i), "r").read()
      temp = f.replace('[', '').replace(']', '').replace(',', '').rstrip().split()
      os.remove(".tmp_%s.dat" %str(start_init + length / np * i))
      #temp = processes[i].communicate()[0].replace('[', '').replace(']', '').replace(',', '').rstrip().split()
      
      for j in range(len(temp) / h_site_count):

        site_occupancy[j + i * len(temp) / h_site_count ] = temp [(j * h_site_count) : ((j + 1) * h_site_count)]


  else: #wird nur ausgefuehrt, fuer Aufruf von subprocess (da np = 1)

    #traj: path to trajectory
    #prmtop: path to parameter file
    #hydration-site: path to hydration-sites as pdb format
    O_coord = MDAnalysis.Universe(h_site).coord
    MDtraj = MDAnalysis.Universe(prmtop, traj)
    density_vs_time_compute (MDtraj, O_coord, site_radius, start_init, stop_init, skip)
    exit()
  data = ''

  for i in site_occupancy.keys():

    for j in site_occupancy[i]:

      data = data + str(j) + ' '

    data = data + '\n'

  o = open('density_vs_time_%s.dat' %name, "w")
  o.write(data)
  o.flush()
  o.close()

  R_input = '''library(bio3d)
pdb_ref <- read.pdb("%s")
ca.inds_ref <- atom.select(pdb_ref, resno=c(1:%s))
data <- read.table("./density_vs_time_%s.dat")
hcor <- cor(data, method="spearman")
view.dccm(hcor, pdb_ref$xyz[ca.inds_ref$xyz], outprefix="%s_corr")
''' %(h_site, str(h_site_count), name, name)
  
  o = open('density_vs_time_%s.R' %name, "w")
  o.write(R_input)
  o.close()

  os.system("R CMD BATCH density_vs_time_%s.R" %name)


  #coord = hydration_site coordinates [[x1,y1,z1], [x2,y2,z2], ...]

def density_vs_time_compute (MDtraj, O_coord, site_radius, start, stop, skip):

  #site_occupancy = {frame#:[occupancy1, occupancy2, ...} The ordering of the occupancies should be exactly the same as in coord.
  site_occupancy = OrderedDict()
  frame_occupancy_init = []
  o = open(".tmp_%s.dat" %str(start), "a")

  for site_count in range(len(O_coord)):
    frame_occupancy_init.append(0)

  for ts in MDtraj.trajectory[start:stop:skip]:
    frame_occupancy = frame_occupancy_init
    for j, site in enumerate(O_coord):
      wat_around_site = MDtraj.selectAtoms("resname WAT and name O and point %s %s %s %s" %(str(site[0]), str(site[1]), str(site[2]), str(site_radius))) 
      frame_occupancy[j] = len(wat_around_site)

    o.write(str(frame_occupancy))
    o.write(', ')
    
    #site_occupancy[ts.frame] = frame_occupancy

  o.close()
  
  return None
  
  #return site_occupancy


def calc_dist (vector1, vector2):

  value = math.sqrt((vector1[0]-vector2[0])**2 + (vector1[1]-vector2[1])**2 + (vector1[2]-vector2[2])**2)
  return value

def calc_mean (vector1, vector2):

  value_vector = []

  value_vector.append((vector1[0] + vector2[0]) / 2.0)
  value_vector.append((vector1[1] + vector2[1]) / 2.0)
  value_vector.append((vector1[2] + vector2[2]) / 2.0)

  return value_vector

def calc_sd (dict_init):

  # dict_init can be real or frac

  scores_str = dict_init.values()
  scores = []

  for i in scores_str:

    if float(i) > 0.0:
      scores.append(float(i))
  sd = np.sqrt(np.var(scores))
  print len(scores)

  if sd == 0.0:
    print "Standard deviation is zero. Check the map."
    exit()

  else:
    print "Standard deviation for all non-zero values is %5.3f." %sd

  dict_new = dict_init

  for key, value in dict_init.items():
    dict_new[key] = str(float(value)/sd)

  return dict_new


def frac2real (dict_init, origin, dim, n, alpha, beta, gamma):

    #dict_init has fractional coordinates!

    dict_real = OrderedDict()
    cos_a = ( math.cos(math.pi*beta/180) * math.cos(math.pi*gamma/180) - math.cos(math.pi*alpha/180) ) / ( math.sin(math.pi*beta/180) *  math.sin(math.pi*gamma/180))
    sin_a = math.sqrt(1.0 - cos_a**2)

    for key, value in dict_init.items():

       X_coord = origin[0] + dim[0] / n[0] * key[0] + dim[1] / n[1] * math.cos(math.pi*gamma/180) * key[1] + dim[2] / n[2] * math.cos(math.pi*beta/180) * key[2]
       Y_coord = origin[1] + 0                      + dim[1] / n[1] * math.sin(math.pi*gamma/180) * key[1] - dim[2] / n[2] * math.sin(math.pi*beta/180) * cos_a * key[2]
       Z_coord = origin[2] + 0                      + 0                                                    + dim[2] / n[2] * math.sin(math.pi*beta/180) * sin_a * key[2]
       dict_real[X_coord, Y_coord, Z_coord] = value

    return dict_real

    
def real_cut (map_data_old, map_data, output):

  bins_old = map_data_old['bins']
  n_old = map_data_old['n']
  dim_old = map_data_old['dim']

  X_bin_size_old = float(dim_old[0]/n_old[0])
  Y_bin_size_old = float(dim_old[1]/n_old[1])
  Z_bin_size_old = float(dim_old[2]/n_old[2])

  old_map_frac = map_data_old['map']

  bins = map_data['bins']
  n = map_data['n']
  dim = map_data['dim']
  origin = map_data['origin']

  alpha = map_data['alpha']
  beta = map_data['beta']
  gamma = map_data['gamma']

  cos_a = ( math.cos(math.pi*beta/180) * math.cos(math.pi*gamma/180) - math.cos(math.pi*alpha/180) ) / ( math.sin(math.pi*beta/180) *  math.sin(math.pi*gamma/180))
  sin_a = math.sqrt(1.0 - cos_a**2)
  crd_real = OrderedDict()

#Find out which values are inside the new box and write them into a new map in real space
  for x in range(0,bins[0]):
    for y in range(0,bins[1]):
      for z in range(0,bins[2]):

        X_crd = origin[0] + dim[0] / n[0] * x + dim[1] / n[1] * math.cos(math.pi*gamma/180) * y + dim[2] / n[2] * math.cos(math.pi*beta/180) * z
        Y_crd = origin[1] + 0                             + dim[1] / n[1] * math.sin(math.pi*gamma/180) * y - dim[2] / n[2] * math.sin(math.pi*beta/180) * cos_a * z
        Z_crd = origin[2] + 0                             + 0                                                               + dim[2] / n[2] * math.sin(math.pi*beta/180) * sin_a * z

        box_checker, basis_vectors = is_in_box([X_crd, Y_crd, Z_crd], map_data_old, True)

        if box_checker:

          X_vector_frac = int(basis_vectors[0] * bins_old[0])
          X_vector_frac_deviation = basis_vectors[0] * float(bins_old[0]) - float(X_vector_frac)

          Y_vector_frac = int(basis_vectors[1] * bins_old[1])
          Y_vector_frac_deviation = basis_vectors[1] * float(bins_old[1]) - float(Y_vector_frac)

          Z_vector_frac = int(basis_vectors[2] * bins_old[2])
          Z_vector_frac_deviation = basis_vectors[2] * float(bins_old[2]) - float(Z_vector_frac)

          if not (0.0 < X_vector_frac_deviation < 1.0 and 0.0 < Y_vector_frac_deviation < 1.0 and 0.0 < Z_vector_frac_deviation < 1.0):

            print "(X Y Z) deviations: %6.5f %6.5f %6.5f" %(X_vector_frac_deviation, Y_vector_frac_deviation, Z_vector_frac_deviation)
            print "(X Y Z) frac-crd used for interpolation: %s %s %s" %(str(X_vector_frac), str(Y_vector_frac), str(Z_vector_frac))
            print "Numerical stuff. Aborting..."
            exit()

          X_value_old_map = float(old_map_frac[X_vector_frac, Y_vector_frac, Z_vector_frac]) * (1 - X_vector_frac_deviation)
          Y_value_old_map = float(old_map_frac[X_vector_frac, Y_vector_frac, Z_vector_frac]) * (1 - Y_vector_frac_deviation)
          Z_value_old_map = float(old_map_frac[X_vector_frac, Y_vector_frac, Z_vector_frac]) * (1 - Z_vector_frac_deviation)

          XX_value_old_map = float(old_map_frac[X_vector_frac+1, Y_vector_frac, Z_vector_frac]) * X_vector_frac_deviation
          YY_value_old_map = float(old_map_frac[X_vector_frac, Y_vector_frac+1, Z_vector_frac]) * Y_vector_frac_deviation
          ZZ_value_old_map = float(old_map_frac[X_vector_frac, Y_vector_frac, Z_vector_frac+1]) * Z_vector_frac_deviation

          crd_real[x, y, z] = (X_value_old_map + Y_value_old_map + Z_value_old_map + XX_value_old_map + YY_value_old_map + ZZ_value_old_map) / 3.0

    print "Finished cycle No %s." %((x+1)*bins[0]**2)

  return crd_real

  
def real2frac (dict_init, origin, dim, n, alpha, beta, gamma):

   nm_dict_frac = OrderedDict()
   cos_a = ( math.cos(math.pi*beta/180) * math.cos(math.pi*gamma/180) - math.cos(math.pi*alpha/180) ) / ( math.sin(math.pi*beta/180) *  math.sin(math.pi*gamma/180))
   sin_a = math.sqrt(1.0 - cos_a**2)

   for x in range(0,n[0]):
    for y in range(0,n[1]):
      for z in range(0,n[2]):

        X_coord = origin[0] + dim[0] / n[0] * x + dim[1] / n[1] * math.cos(math.pi*gamma/180) * y + dim[2] / n[2] * math.cos(math.pi*beta/180) * z
        Y_coord = origin[1] + 0                 + dim[1] / n[1] * math.sin(math.pi*gamma/180) * y - dim[2] / n[2] * math.sin(math.pi*beta/180) * cos_a * z
        Z_coord = origin[2] + 0                 + 0                                               + dim[2] / n[2] * math.sin(math.pi*beta/180) * sin_a * z
        nm_dict_frac[x,y,z] = 'NA'

        for key, value in dict_init.items():

          if abs(calc_dist([X_coord, Y_coord, Z_coord], [key[0], key[1], key[2]])) < abs(calc_dist([0.0,0.0,0.0], [dim[0] / n[0], dim[0] / n[0], dim[0] / n[0]]))/2:

            nm_dict_frac[x,y,z] = str(value)

   return nm_dict_frac

def write_dx (dict_init, origin, dim, new_map, end_bins, start_bins=[0,0,0], factor=1.0):

  #dict_init must be frac_coordinates!

  bins = [int(end_bins[0])-start_bins[0], int(end_bins[1])-start_bins[1], int(end_bins[2])-start_bins[2]]

  header = '''object 1 class gridpositions counts %s %s %s
origin %s %s %s
delta %s 0 0
delta 0 %s 0
delta 0 0 %s
object 2 class gridconnections counts %s %s %s
object 3 class array type float rank 0 items %s data follows
''' %(str(bins[0]), str(bins[1]), str(bins[2]), str(origin[0]), str(origin[1]), str(origin[2]), str(dim[0]/bins[0]), str(dim[1]/bins[1]), str(dim[2]/bins[2]), str(bins[0]), str(bins[1]), str(bins[2]), str(bins[2] * bins[1] * bins[0]))

  data = ''

  i = 0

  for x_frac in range(start_bins[0], int(end_bins[0])):
    for y_frac in range(start_bins[1], int(end_bins[1])):
      for z_frac in range(start_bins[2], int(end_bins[2])):

        data = data + str(float(dict_init[x_frac, y_frac, z_frac]) * factor) + ' '
        i = i + 1

        if i == 3:

           data = data + '\n'
           i = 0

  o = open('%s.dx' %new_map, "w")
  o.write(header + data)
  o.close()


def map_diff (dict_init_1, dict_init_2, bins1, factor1=1.0, factor2=-1.0):
#Both maps dict_init_1 and dict_init_2 must be in fractional coordinates

   dict_diff = OrderedDict()
   for x_frac in range (0, bins1[0]):
     for y_frac in range (0, bins1[1]):
       for z_frac in range (0, bins1[2]):
         dict_diff[x_frac, y_frac, z_frac] = factor1*float(dict_init_1[x_frac, y_frac, z_frac]) + factor2*float(dict_init_2[x_frac, y_frac, z_frac])
   return dict_diff


def refine (dict_init, max_dist, Nr_nearest_neighbours=0, map_data='empty'):

  #This method needs dict_init in real coordinates. Map_data is the "raw" initial data with fractional coordinates.

  if map_data != 'empty' and Nr_nearest_neighbours > 0:

    bins = map_data['bins']
    n = map_data['n']
    dim = map_data['dim']

  weighted_best_of = OrderedDict()

  for crd_ref in dict_init.keys():

     best_of = OrderedDict()

     if Nr_nearest_neighbours > 0:

       for i in range(0, Nr_nearest_neighbours):

#         shortest_distance = []

#         shortest_distance.append(dim[0])
#         shortest_distance.append(dim[1])
#         shortest_distance.append(dim[2])

         shortest_distance = math.sqrt( 3 * dim[0]**2 )

         for crd in dict_init.keys():

           if crd not in best_of:

#            dist = []

#            dist.append(abs(crd[0] - crd_ref[0]))
#            dist.append(abs(crd[1] - crd_ref[1]))
#            dist.append(abs(crd[2] - crd_ref[2]))

#            if calc_dist([dist[0], dist[1], dist[2]],[0.0, 0.0, 0.0]) < calc_dist([shortest_distance[0], shortest_distance[1], shortest_distance[2]], [0.0, 0.0, 0.0]) and dist[0] <= max_dist[0] and dist[1] <= max_dist[1] and dist[2] <= max_dist[2]:

#              shortest_distance = dist
#              best_of[crd] = dict_init[crd]

            dist = math.sqrt( (crd[0]-crd_ref[0])**2 + (crd[1]-crd_ref[1])**2 + (crd[2]-crd_ref[2])**2 )

            if dist <= shortest_distance and dist < max_dist:

              shortest_distance = dist
              best_of[crd] = dict_init[crd]

     else:

       for crd in dict_init.keys():

#         dist = []

#         dist.append(abs(crd[0] - crd_ref[0]))
#         dist.append(abs(crd[1] - crd_ref[1]))
#         dist.append(abs(crd[2] - crd_ref[2]))

#         if crd not in best_of and dist[0] <= max_dist[0] and dist[1] <= max_dist[1] and dist[2] <= max_dist[2]:

#           best_of[crd] = dict_init[crd]

           if crd not in best_of:

             dist = math.sqrt( (crd[0]-crd_ref[0])**2 + (crd[1]-crd_ref[1])**2 + (crd[2]-crd_ref[2])**2 )

             if dist <= max_dist:

               best_of[crd] = dict_init[crd]

     coord_average = []
     coord_average.append(0.0)
     coord_average.append(0.0)
     coord_average.append(0.0)

     best_of_valuesum = 0.0

     if len(best_of) == 0: best_of[crd_ref] = dict_init[crd_ref]

     for key in best_of.keys():

        coord_average[0] = coord_average[0] + key[0] * float(best_of[key])
        coord_average[1] = coord_average[1] + key[1] * float(best_of[key])
        coord_average[2] = coord_average[2] + key[2] * float(best_of[key])
        best_of_valuesum = best_of_valuesum + float(best_of[key])

     coord_average[0] = coord_average[0]/best_of_valuesum
     coord_average[1] = coord_average[1]/best_of_valuesum
     coord_average[2] = coord_average[2]/best_of_valuesum

     weighted_best_of[coord_average[0], coord_average[1], coord_average[2]] = best_of_valuesum / float(len(best_of))

  return weighted_best_of



def find_maximum (map_data, minimum_density, refine_factor=1.0):

   bins = map_data['bins']
   dict_init = map_data['map']
   n = map_data['n']
   dim = map_data['dim']

   maximum = OrderedDict()

   for x in range(0, bins[0]):
     for y in range(0, bins[1]):
       for z in range(0, bins[2]):
         if float(dict_init[x,y,z]) > float(minimum_density):
           maximum[x,y,z] = dict_init[x,y,z]
   maximum = frac2real(maximum, map_data['origin'], map_data['dim'], map_data['n'], map_data['alpha'], map_data['beta'], map_data['gamma'])

   #write_pdb(maximum, 'maximum')

#The max_dist_factor is very important!  An increase will lead to averaging of hydration sites that are to far away from each other, a decrease will lead to too many hydration sites in direct neighborhood.

#   max_dist_factor = []

#   max_dist_factor.append(float(dim[0] / n[0]))
#   max_dist_factor.append(float(dim[1] / n[1]))
#   max_dist_factor.append(float(dim[2] / n[2]))

   print "Initial preselection of six nearest neighbors within delta = %s A." %refine_factor
   refine_map1 = refine (maximum, refine_factor, 6, map_data)

   #write_pdb(refine_map1, 'refine_best_of_6')

#This Method ensures that all hydration sites are minimum 1.0 A apparte from each other
#The results are better compared to the other method (see below)

#   max_dist_factor = []

#   max_dist_factor.append(float(dim[0] / n[0] * refine_factor))
#   max_dist_factor.append(float(dim[1] / n[1] * refine_factor))
#   max_dist_factor.append(float(dim[2] / n[2] * refine_factor))

   refine_map2 = refine_map1
   refine_map2_old = OrderedDict()
   first = True

   while (len(refine_map2) < len(refine_map2_old) or first):

     input1 = len(refine_map2)
     refine_map2_old = refine_map2
     refine_map2 = refine (refine_map2, refine_factor, 0)
     first = False
     input2 = len(refine_map2)
     print 'Refinement of %s grid points' %str(input1 - input2)

   return refine_map2

# This method practically removes "bad contacts" of the second refined map 'refine_map2'
# It cannot remove all bad contacts, so the other method is better.

#   for i in range(1,10):
#
#    refine_map3 = refine (refine_map2, (1.0 - i * 0.1), 0, map_data, minimum_density)
#
#    refine_map2 = refine_map3
#
#   return refine_map3

def find_max_dc (dict_init, hradius, threshold):

  crds   = []
  values = []

#build numpy version of map
  for i, crd in enumerate(dict_init.keys()):
#filter
    crds.append(crd)
    values.append(dict_init[crd])

  crds     = np.array(crds)
  values   = np.array(values)

  out_crds = {}

  def update_crds(p, v, crds, values):
    out_crds[p[0], p[1], p[2]] = v

    updater = [ np.where( distance.cdist(p[np.newaxis,:], crds, 'euclidean') > hradius )[1] ]
    crds    = crds[updater]
    values  = values[updater]

    return crds, values
  
  while len(np.where (values > threshold)[0]) > 0:

    new = np.where ( values == np.amax(values) )[0]
    crds, values = update_crds(crds[new][0], values[new][0], crds, values)

  return out_crds


def find_max_dc_old (dict_init, hradius, threshold):

  #dict_init must be real

  dict_max = OrderedDict()
  dict_ = dict_init
  while True:

    dict_old = dict_
    dict_ = OrderedDict()
    max_value = 0.0
    max_key = [0.0, 0.0, 0.0]

    for key, value in dict_old.items():

      if float(value) > max_value and float(value) > threshold: 

        max_value = float(value)

        max_key[0] = key[0]
        max_key[1] = key[1]
        max_key[2] = key[2]

    if max_value == 0.0:

      print "Found %s maxima." %len(dict_max)
      return dict_max

    else:

      dict_max[max_key[0], max_key[1], max_key[2]] = str(max_value)

    for key, value in dict_old.items():

      if calc_dist(key, max_key) < hradius:

        dict_[key] = '0.0'

      else:

        dict_[key] = value


def write_pdb (dict_init, new_map, dict_indexing={}):

  #dict_init must be real coordinates!

  data = 'REMARK default_name\n'

  if dict_indexing == {}:

    for i, key in enumerate(dict_init.keys()):
      data = data + '%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' %('HETATM',i+1,'O','', 'HOH', 'A', i+1, '', key[0], key[1], key[2], 0.00, float(dict_init[key]), 'O', '')

  else:

    for key in dict_init.keys():
      data = data + '%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' %('HETATM',dict_indexing[key],'O','', 'HOH', 'A', dict_indexing[key], '', key[0], key[1], key[2], 0.00, float(dict_init[key]), 'O', '')


  data = data + 'END\n'

  o = open('%s.pdb' %new_map, "w")
  o.write(data)
  o.close()


def read_pdb (file):

  crd_file_ref = open(file,"r")
  crd_file = crd_file_ref.readlines()
  crd_file_ref.close()

  dict_init = OrderedDict()

  key = [0.0, 0.0, 0.0]

  for i, item in enumerate(crd_file):

    if not (item[0:6].rstrip() == 'ATOM' or item[0:6].rstrip() == 'HETATM'):
      continue

    if i <= 9999:

      key[0] = float(item.rstrip()[30:38])
      key[1] = float(item.rstrip()[38:46])
      key[2] = float(item.rstrip()[46:54])

# load B-Factors as values
      dict_init[key[0], key[1], key[2]] = item.rstrip()[54:59]

    if 9999 < i <= 99999:

      key[0] = float(item.rstrip()[31:39])
      key[1] = float(item.rstrip()[39:47])
      key[2] = float(item.rstrip()[47:55])

# load B-Factors as values
      dict_init[key[0], key[1], key[2]] = item.rstrip()[55:60]

    if i > 99999:

      key[0] = float(item.rstrip()[33:41])
      key[1] = float(item.rstrip()[41:49])
      key[2] = float(item.rstrip()[49:57])

# load B-Factors as values
      dict_init[key[0], key[1], key[2]] = item.rstrip()[57:62]

  return dict_init


def write_map_table(map_data, output_name, additional_data={}):

  dict_map = map_data['map']
  bins = map_data['bins']

  length_additional_data = len(additional_data)

  if length_additional_data != 0 and length_additional_data != len(dict_map.values()):

    print "Something wrong with map-data and additional data."
    exit()

  if length_additional_data == len(dict_map.values()):

    data = 'x_crd y_crd z_crd value'

    for i in range(0,length_additional_data):
      data = data + ' add_value%d' %i

    data = data + '\n'

    for x in range(0,bins[0]):
      for y in range(0,bins[1]):
        for z in range(0,bins[2]):

          data = data + '%d %d %d %s %s\n' %(x,y,z, dict_map[x,y,z], str(additional_data[x,y,z]))

  if length_additional_data == 0:

    data = 'x_crd y_crd z_crd value\n'

    for x in range(0,bins[0]):
      for y in range(0,bins[1]):
        for z in range(0,bins[2]):

          data = data + '%d %d %d %s\n' %(x,y,z, dict_map[x,y,z],)

  o = open('%s.dat' %output_name, "w")
  o.write(data)
  o.close()



def write_table (dict_map, output_name, additional_data={}):

  data = 'x_crd y_crd z_crd value\n'

  length_additional_data = len(additional_data)

  if length_additional_data != 0 and length_additional_data != len(dict_map.values()):

    print "Something wrong with map-data and additional data."
    exit()

  if length_additional_data == len(dict_map.values()):

    for key, value in dict_map.items():

      data = data + '%s %s %s %s %s\n' %(str(key[0]), str(key[1]), str(key[2]), value, str(additional_data[key[0], key[1], key[2]]))

  if length_additional_data == 0:

    for key, value in dict_map.items():

      data = data + '%s %s %s %s\n' %(str(key[0]), str(key[1]), str(key[2]), value)

  o = open('%s.dat' %output_name, "w")
  o.write(data)
  o.close()


def calc_radial_density_distribution (map_data, reference_model, output_name, reference_model_selection="all", threshold=0.0):

  alpha = map_data['alpha']
  beta = map_data['beta']
  gamma = map_data['gamma']

  origin = map_data['origin']
  bins = map_data['bins']
  n = map_data['n']
  dim = map_data['dim']

  np_map = init2np(map_data['map'], bins)

  new_dict = OrderedDict()

# selection like "resname WAT and name O"

  reference_model = MDAnalysis.Universe(reference_model).selectAtoms(reference_model_selection).atoms.coordinates()

#x,y,z coordinates via reference_model[0],[1],[2]

  raster_size = dim[2] / n[2]
  raster_length = int(math.sqrt(bins[0]**2 + bins[1]**2 + bins[2]**2))

  if dim[1] / n[1] > raster_size:

    raster_size = dim[1] / n[1]

  if dim[0] / n[0] > raster_size:

    raster_size = dim[0] / n[0]

  dist_vs_dens = []

  for i in range(0, raster_length): dist_vs_dens.append(0.0)

# value_feeding

  cos_a = ( math.cos(math.pi*beta/180) * math.cos(math.pi*gamma/180) - math.cos(math.pi*alpha/180) ) / ( math.sin(math.pi*beta/180) *  math.sin(math.pi*gamma/180))
  sin_a = math.sqrt(1.0 - cos_a**2)

  for x in range(0, bins[0]):
    for y in range(0, bins[1]):
      for z in range(0, bins[2]):

        X_coord = origin[0] + dim[0] / n[0] * x + dim[1] / n[1] * math.cos(math.pi*gamma/180) * y + dim[2] / n[2] * math.cos(math.pi*beta/180) * z
        Y_coord = origin[1] + 0                 + dim[1] / n[1] * math.sin(math.pi*gamma/180) * y - dim[2] / n[2] * math.sin(math.pi*beta/180) * cos_a * z
        Z_coord = origin[2] + 0                 + 0                                               + dim[2] / n[2] * math.sin(math.pi*beta/180) * sin_a * z

        for site in reference_model:

          dist = calc_dist([site[0], site[1], site[2]], [X_coord, Y_coord, Z_coord])

          if is_in_box(site, map_data) and np_map[x][y][z] > threshold:

            dist_vs_dens[int(dist / raster_size)] = dist_vs_dens[int(dist / raster_size)] + np_map[x][y][z] * 0.0334 / (4.0 * np.pi * dist**2 * raster_size)

    print "Finished cycle No %s." %((x+1)*bins[0]**2)

#transfer values to readable-format

  data = 'distance mean_value\n'

  for i in range(0, raster_length):

    data = data + '%9.3f %10.5f' %(i*raster_size, dist_vs_dens[i])
    data = data + '\n'

  o = open('%s.dat' %output_name, "w")
  o.write(data)
  o.close()


def is_in_box(crd, map_data,return_coeffs=False):

  alpha = map_data['alpha']
  beta = map_data['beta']
  gamma = map_data['gamma']

  origin = map_data['origin']
  bins = map_data['bins']
  dim = map_data['dim']
  n = map_data['n']

  check = False

  cos_a = ( math.cos(math.pi*beta/180) * math.cos(math.pi*gamma/180) - math.cos(math.pi*alpha/180) ) / ( math.sin(math.pi*beta/180) * math.sin(math.pi*gamma/180))

  sin_a = math.sqrt(1.0 - cos_a**2)

  x_max = []
  y_max = []
  z_max = []

  x_max.append(dim[0] / n[0] * bins[0] + dim[1] / n[1] * math.cos(math.pi*gamma/180) * 0 + dim[2] / n[2] * math.cos(math.pi*beta/180) * 0)
  x_max.append(0                       + dim[1] / n[1] * math.sin(math.pi*gamma/180) * 0 - dim[2] / n[2] * math.sin(math.pi*beta/180) * cos_a * 0)
  x_max.append(0                       + 0                                                     + dim[2] / n[2] * math.sin(math.pi*beta/180) * sin_a * 0)

  y_max.append(dim[0] / n[0] * 0 + dim[1] / n[1] * math.cos(math.pi*gamma/180) * bins[1] + dim[2] / n[2] * math.cos(math.pi*beta/180) * 0)
  y_max.append(0                       + dim[1] / n[1] * math.sin(math.pi*gamma/180) * bins[1] - dim[2] / n[2] * math.sin(math.pi*beta/180) * cos_a * 0)
  y_max.append(0                       + 0                                                     + dim[2] / n[2] * math.sin(math.pi*beta/180) * sin_a * 0)

  z_max.append(dim[0] / n[0] * 0 + dim[1] / n[1] * math.cos(math.pi*gamma/180) * 0 + dim[2] / n[2] * math.cos(math.pi*beta/180) * bins[2])
  z_max.append(0                       + dim[1] / n[1] * math.sin(math.pi*gamma/180) * 0 - dim[2] / n[2] * math.sin(math.pi*beta/180) * cos_a * bins[2])
  z_max.append(0                       + 0                                                     + dim[2] / n[2] * math.sin(math.pi*beta/180) * sin_a * bins[2])

#Solve ''crd[0] = origin[0] + a * (x_max) + b * (- origin[1]) + c * (- origin[2]), crd[1] = ...''
  A = np.array([[x_max[0],y_max[0],z_max[0]], [x_max[1],y_max[1],z_max[1]], [x_max[2],y_max[2],z_max[2]]])
  b = np.array([crd[0]-origin[0], crd[1]-origin[1], crd[2]-origin[2]])
  solve = np.linalg.solve(A, b)

  if 0.0 <= solve[0] <= 1.0 and 0.0 <= solve[1] <= 1.0 and 0.0 <= solve[2] <= 1.0:

    check = True

  if return_coeffs == False:

    return check

  else:

    return check, solve


def site_site (init_map, output_name):

#returns under triangle of water-water-distance matrix
  new_map = {}
  data = ''
  length = len(init_map.keys())

  for i, key1 in enumerate(init_map.keys()):
    for j in range(i):
      data = data + 'NA '
    for key2 in init_map.keys()[(i+1):length]:
      dist = calc_dist(key1, key2)
      data = data + '%5.3f' %dist + ' '
    data = data + '\n'

  o = open('%s.dat' %output_name, "w")
  o.write(data)
  o.close()


def nearest_neighbour_RMSD (sites, ref_sites, map_data):

#both site maps must be real.

  RMSD = []

  for ref_site in ref_sites.keys():

    nearest_neighbour_distance = 100000000.00

    if is_in_box([ref_site[0], ref_site[1], ref_site[2]], map_data):

      for site in sites.keys():

        dist = calc_dist([site[0],site[1],site[2]], [ref_site[0], ref_site[1], ref_site[2]])

        if dist < nearest_neighbour_distance:

          nearest_neighbour_distance = dist

      if nearest_neighbour_distance < 5.0:

        RMSD.append(nearest_neighbour_distance**2)

  RMSD_value = 0.0

  for i in RMSD:

    RMSD_value = RMSD_value + i

  return 1.0/len(RMSD) * math.sqrt(RMSD_value)

def nearest_neighbour_dist (sites, ref_sites, map_data, threshold):

  site_data = {}

#Default value declaration...
  if threshold==None: threshold=5.0
  for ref_site in ref_sites.keys():

    nearest_neighbour_distance = 100000000.00

    if is_in_box([ref_site[0], ref_site[1], ref_site[2]], map_data):

      for site in sites.keys():

        dist = calc_dist([site[0],site[1],site[2]], [ref_site[0], ref_site[1], ref_site[2]])

        if dist < nearest_neighbour_distance:

          nearest_neighbour_distance = dist

      site_data[ref_site] = nearest_neighbour_distance

  return site_data


def find_max_grad (map_data):

#fractional map
  init_map = map_data['map']

  alpha = map_data['alpha']
  beta = map_data['beta']
  gamma = map_data['gamma']

  origin = map_data['origin']
  bins = map_data['bins']
  n = map_data['n']
  dim = map_data['dim']

  np_map = init2np(init_map, bins)

  new_init_map = np2init(np_map, bins)

  x_grad, y_grad, z_grad = np.gradient(np.array(np_map, dtype=np.float), edge_order=2)

  x_grad_init = np2init(x_grad, bins)
  y_grad_init = np2init(y_grad, bins)
  z_grad_init = np2init(z_grad, bins)

#Calc and Feed it with values

  new_map = {}

  cos_a = ( math.cos(math.pi*beta/180) * math.cos(math.pi*gamma/180) - math.cos(math.pi*alpha/180) ) / ( math.sin(math.pi*beta/180) *  math.sin(math.pi*gamma/180))

  sin_a = math.sqrt(1.0 - cos_a**2)

  for x in range(1, bins[0]-1):
    for y in range(1, bins[1]-1):
      for z in range(1, bins[2]-1):

        differential = x_grad[x][y][z] + y_grad[x][y][z] + z_grad[x][y][z]

        if 0.0 < x_grad[x-1][y][z] and 0.0 > x_grad[x+1][y][z] and 0.0 < y_grad[x][y-1][z] and 0.0 > y_grad[x][y+1][z] and 0.0 < z_grad[x][y][z-1] and 0.0 > z_grad[x][y][z+1]:

          X_coord = origin[0] + dim[0] / n[0] * x + dim[1] / n[1] * math.cos(math.pi*gamma/180) * y + dim[2] / n[2] * math.cos(math.pi*beta/180) * z
          Y_coord = origin[1] + 0                 + dim[1] / n[1] * math.sin(math.pi*gamma/180) * y - dim[2] / n[2] * math.sin(math.pi*beta/180) * cos_a * z
          Z_coord = origin[2] + 0                 + 0                                               + dim[2] / n[2] * math.sin(math.pi*beta/180) * sin_a * z
          new_map[X_coord, Y_coord, Z_coord] = init_map[x,y,z]

  print "Found %s values." %len(new_map)

  return new_map


def init2np(init_map, bins):

  x = np.array(range(0,bins[0]), dtype=np.float)
  y = np.array(range(0,bins[1]), dtype=np.float)
  z = np.array(range(0,bins[2]), dtype=np.float)

  xx, yy, zz = np.meshgrid(y, x, z)

  value = xx

  for x in range(0,bins[0]):
    for y in range(0,bins[1]):
      for z in range(0,bins[2]):

        value[x][y][z] = float(init_map[x,y,z])

  return value

def np2init (np_map, bins):

  init_map = {}

  for x in range(0,bins[0]):
    for y in range(0,bins[1]):
      for z in range(0,bins[2]):

        init_map[x,y,z] = str(np_map[x][y][z])

  return init_map

def smooth_gauss(init_map, bins):

  np_map = init2np(init_map, bins)
  np_map = ndimage.gaussian_filter(np_map,1)
  new_map = np2init(np_map, bins)

  return new_map

def find_max_gradient_preselection(map_data, refine_factor=1.0):

  bins = map_data['bins']
  dict_init = map_data['map']
  n = map_data['n']
  dim = map_data['dim']

  max_grad_map = find_max_grad(map_data)

  max_dist_factor = []

  max_dist_factor.append(float(dim[0] / n[0] * refine_factor))
  max_dist_factor.append(float(dim[1] / n[1] * refine_factor))
  max_dist_factor.append(float(dim[2] / n[2] * refine_factor))

  refine_map2 = max_grad_map

  refine_map2_old = {}

  first = True

  while (len(refine_map2) < len(refine_map2_old) or first):

    input1 = len(refine_map2)
    refine_map2_old = refine_map2
    refine_map2 = refine (refine_map2, max_dist_factor, 0, map_data)
    first = False
    input2 = len(refine_map2)
    print 'Refinement of %s grid points' %str(input1 - input2)

  return refine_map2



def placevent(map_data, new_map):

  #dict_init must be frac and a population distribution, no g_O !

  dict_old = map_data['map']
  dim = map_data['dim']
  n = map_data['n']
  bins = map_data['bins']

  dict_max = {}

  B_factor_map = {}

  while True:

    density_round_maximum = 0.0
    max_key = [0, 0, 0]


# first find the maximum of the actual population distribution

    for key, value in dict_old.items():

      if float(value) > density_round_maximum: 

        density_round_maximum = float(value)
        max_key[0] = key[0]
        max_key[1] = key[1]
        max_key[2] = key[2]

    if density_round_maximum == 0.0:

      print "Found %s maxima." %len(dict_max)
      real_map = frac2real (dict_max, map_data['origin'], dim, n, map_data['alpha'], map_data['beta'], map_data['gamma'])
      data = 'REMARK default_name\n'

      for i, key in enumerate(real_map.keys()):

        data = data + '%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\nTER\n' %('HETATM',i+1,'MA','', 'MAP', 'A', 1, '', key[0], key[1], key[2], float(real_map[key][0]), float(real_map[key][1]), 'O', '')

      data = data + 'END\n'

      o = open('%s.pdb' %new_map, "w")
      o.write(data)
      o.close()

      return

#Now lets cut out all values in radius around maximum, such that the total population in this radius equals one

    i = 1

    while density_round_maximum <= 1.0 and 0 < max_key[0]+i < bins[0] and 0 < max_key[0]-i < bins[0] and 0 < max_key[1]+i < bins[1] and 0 < max_key[1]-i < bins[1] and 0 < max_key[2]+i < bins[2] and 0 < max_key[2]-i < bins[2]:

#all nearest neighbors...      
      density_round_maximum = density_round_maximum + float(dict_old[max_key[0]+i, max_key[1] , max_key[2]]) + float(dict_old[max_key[0]-i, max_key[1], max_key[2]]) + float(dict_old[max_key[0], max_key[1]+i, max_key[2]]) + float(dict_old[max_key[0], max_key[1]-i, max_key[2]]) + float(dict_old[max_key[0], max_key[1] , max_key[2]+i]) + float(dict_old[max_key[0], max_key[1], max_key[2]-i])
      dict_old[max_key[0]+i, max_key[1], max_key[2]] = '0.0'
      dict_old[max_key[0]-i, max_key[1], max_key[2]] = '0.0'
      dict_old[max_key[0], max_key[1]+i, max_key[2]] = '0.0'
      dict_old[max_key[0], max_key[1]-i, max_key[2]] = '0.0'
      dict_old[max_key[0], max_key[1], max_key[2]+i] = '0.0'
      dict_old[max_key[0], max_key[1], max_key[2]-i] = '0.0'
      smallest_radius_round_maximum = (i*dim[0]/n[0] + i*dim[1]/n[1] + i*dim[2]/n[2]) / 3.0

#and now all plane diagonal neighbors
      if (density_round_maximum <= 1.0):
        density_round_maximum = density_round_maximum + float(dict_old[max_key[0]+i, max_key[1]+i , max_key[2]]) + float(dict_old[max_key[0]+i, max_key[1]-i , max_key[2]]) + float(dict_old[max_key[0]-i, max_key[1]+i, max_key[2]]) + float(dict_old[max_key[0]-i, max_key[1]-i, max_key[2]]) + float(dict_old[max_key[0]+i, max_key[1], max_key[2]+i]) + float(dict_old[max_key[0]+i, max_key[1], max_key[2]-i]) + float(dict_old[max_key[0]-i, max_key[1], max_key[2]+i]) + float(dict_old[max_key[0]-i, max_key[1], max_key[2]-i]) + float(dict_old[max_key[0], max_key[1]+i, max_key[2]+i]) + float(dict_old[max_key[0], max_key[1]+i, max_key[2]-i]) + float(dict_old[max_key[0], max_key[1]-i, max_key[2]+i]) + float(dict_old[max_key[0], max_key[1]-i, max_key[2]-i])
        dict_old[max_key[0]+i, max_key[1]+i, max_key[2]] = '0.0'
        dict_old[max_key[0]+i, max_key[1]-i, max_key[2]] = '0.0'
        dict_old[max_key[0]-i, max_key[1]+i, max_key[2]] = '0.0'
        dict_old[max_key[0]-i, max_key[1]-i, max_key[2]] = '0.0'
        dict_old[max_key[0]+i, max_key[1], max_key[2]+i] = '0.0'
        dict_old[max_key[0]+i, max_key[1], max_key[2]-i] = '0.0'
        dict_old[max_key[0]-i, max_key[1], max_key[2]+i] = '0.0'
        dict_old[max_key[0]-i, max_key[1], max_key[2]-i] = '0.0'
        dict_old[max_key[0], max_key[1]+i, max_key[2]+i] = '0.0'
        dict_old[max_key[0], max_key[1]+i, max_key[2]-i] = '0.0'
        dict_old[max_key[0], max_key[1]-i, max_key[2]+i] = '0.0'
        dict_old[max_key[0], max_key[1]-i, max_key[2]-i] = '0.0'
        smallest_radius_round_maximum = (calc_dist([i*dim[0]/n[0], i*dim[1]/n[1], 0.0], [0.0, 0.0, 0.0]) + calc_dist([i*dim[0]/n[0], 0.0, i*dim[2]/n[2]], [0.0, 0.0, 0.0]) + calc_dist([0.0, i*dim[1]/n[1], i*dim[2]/n[2]], [0.0, 0.0, 0.0]))


#and now all cubic diagonal neighbors
      if (density_round_maximum <= 1.0):
        density_round_maximum = density_round_maximum + float(dict_old[max_key[0]+i, max_key[1]+i, max_key[2]+i]) + float(dict_old[max_key[0]+i, max_key[1]+i , max_key[2]-i]) + float(dict_old[max_key[0]+i, max_key[1]-i , max_key[2]-i]) + float(dict_old[max_key[0]+i, max_key[1]-i , max_key[2]+i]) + float(dict_old[max_key[0]-i, max_key[1]+i , max_key[2]+i]) + float(dict_old[max_key[0]-i, max_key[1]+i , max_key[2]-i]) + float(dict_old[max_key[0]-i, max_key[1]-i , max_key[2]-i]) + float(dict_old[max_key[0]-i, max_key[1]-i , max_key[2]+i])
        dict_old[max_key[0]+i, max_key[1]+i , max_key[2]+i] = '0.0' 
        dict_old[max_key[0]+i, max_key[1]+i , max_key[2]-i] = '0.0'
        dict_old[max_key[0]+i, max_key[1]-i , max_key[2]-i] = '0.0'
        dict_old[max_key[0]+i, max_key[1]-i , max_key[2]+i] = '0.0'
        dict_old[max_key[0]-i, max_key[1]+i , max_key[2]+i] = '0.0'
        dict_old[max_key[0]-i, max_key[1]+i , max_key[2]-i] = '0.0'
        dict_old[max_key[0]-i, max_key[1]-i , max_key[2]-i] = '0.0'
        dict_old[max_key[0]-i, max_key[1]-i , max_key[2]+i] = '0.0'
        smallest_radius_round_maximum = calc_dist([dim[0]/n[0], dim[1]/n[1], dim[2]/n[2]], [0.0, 0.0, 0.0])

      i = i + 1

    dict_max[max_key[0], max_key[1], max_key[2]] = (str(density_round_maximum), str(8.0 * math.pi**2 * smallest_radius_round_maximum))
    dict_old[max_key[0], max_key[1], max_key[2]] = '0.0'



def find_random(map_data, no_of_sites):

#fractional coordinates

  new_crd = {}

  bins = map_data['bins']
  n = map_data['n']
  dim = map_data['dim']
  origin = map_data['origin']

  dict_init = map_data['map']

  alpha = map_data['alpha']
  beta = map_data['beta']
  gamma = map_data['gamma']

  cos_a = ( math.cos(math.pi*beta/180) * math.cos(math.pi*gamma/180) - math.cos(math.pi*alpha/180) ) / ( math.sin(math.pi*beta/180) *  math.sin(math.pi*gamma/180))
  sin_a = math.sqrt(1.0 - cos_a**2)
  i = 0
  
  while i < no_of_sites:

    x_random = random.random() * (bins[0] - 1)
    y_random = random.random() * (bins[1] - 1)
    z_random = random.random() * (bins[2] - 1)
    
    X_crd = origin[0] + dim[0] / n[0] * x_random + dim[1] / n[1] * math.cos(math.pi*gamma/180) * y_random + dim[2] / n[2] * math.cos(math.pi*beta/180) * z_random
    Y_crd = origin[1] + 0                        + dim[1] / n[1] * math.sin(math.pi*gamma/180) * y_random - dim[2] / n[2] * math.sin(math.pi*beta/180) * cos_a * z_random
    Z_crd = origin[2] + 0                        + 0                                                      + dim[2] / n[2] * math.sin(math.pi*beta/180) * sin_a * z_random

    if float(dict_init[int(round(x_random, 0)), int(round(y_random, 0)), int(round(z_random, 0))]) > 1.00:

      new_crd[X_crd, Y_crd, Z_crd] = '0.0'
      i = i + 1

  return new_crd

def boltzmann_inverter(init_dict_frac, bins):

  dict_inverted = {}

  for x in range(0, bins[0]):
    for y in range(0, bins[1]):
      for z in range(0, bins[2]):

        dict_inverted[x,y,z] = np.exp(-float(init_dict_frac[x,y,z]) / (0.0019872041 * 300)) # Gaskonstante(kcal/mol*K) * T

  return dict_inverted

def who_is_in_the_box(pdb_crd, box):

  is_in_box_string = ''

  for i, crd in enumerate(pdb_crd.keys()):

    if is_in_box(crd, box):

      is_in_box_string = is_in_box_string + "%d," %i

  print is_in_box_string


def map_gist(gist_data, Eww_matrix, reference_energy,  calculation_density, reference_density, output):

  Eww_norm       = gist_data['Eww_norm_unref']

  gO             = gist_data['gO']

  bins           = gist_data['bins']
  n              = gist_data['n']
  dim            = gist_data['dim']

  voxel_volume   = (dim[0] / n[0]) * (dim[1] / n[1]) * (dim[2] / n[2])
  density = {}

  for x in range(0, bins[0]):
    for y in range(0, bins[1]):
      for z in range(0, bins[2]):

        density[x,y,z] = gO[x,y,z] * reference_density

        if (x,y,z) in Eww_matrix:

          Eww_norm[x,y,z] = 2 * Eww_norm[x,y,z] - (2 * sum(Eww_matrix[x,y,z]) / (voxel_volume * density[x,y,z]))

  dH  = map_diff(Eww_norm,                        gist_data['Esw_norm'],      gist_data['bins'], 1.0,  1.0)
  TdS = map_diff(gist_data['dTSorient_norm'],     gist_data['dTStrans_norm'], gist_data['bins'], 1.0,  1.0)
  dG  = map_diff(dH,                              TdS,                        gist_data['bins'], 1.0, -1.0)
  
  write_dx(dH,                          gist_data['origin'], gist_data['dim'], '%s_dH'             %output, gist_data['bins'])
  write_dx(TdS,                         gist_data['origin'], gist_data['dim'], '%s_TdS'            %output, gist_data['bins'])
  write_dx(dG,                          gist_data['origin'], gist_data['dim'], '%s_dG'             %output, gist_data['bins'])
  write_dx(density,                     gist_data['origin'], gist_data['dim'], '%s_density'        %output, gist_data['bins'])
  write_dx(gist_data['pop'],            gist_data['origin'], gist_data['dim'], '%s_pop'            %output, gist_data['bins'])
  write_dx(gist_data['gO'],             gist_data['origin'], gist_data['dim'], '%s_gO'             %output, gist_data['bins'])
  write_dx(gist_data['gH'],             gist_data['origin'], gist_data['dim'], '%s_gH'             %output, gist_data['bins'])
  write_dx(gist_data['dTStrans_dens'],  gist_data['origin'], gist_data['dim'], '%s_dTStrans_dens'  %output, gist_data['bins'])
  write_dx(gist_data['dTStrans_norm'],  gist_data['origin'], gist_data['dim'], '%s_dTStrans_norm'  %output, gist_data['bins'])
  write_dx(gist_data['dTSorient_dens'], gist_data['origin'], gist_data['dim'], '%s_dTSorient_dens' %output, gist_data['bins'])
  write_dx(gist_data['dTSorient_norm'], gist_data['origin'], gist_data['dim'], '%s_dTSorient_norm' %output, gist_data['bins'])
  write_dx(gist_data['Esw_dens'],       gist_data['origin'], gist_data['dim'], '%s_Esw_dens'       %output, gist_data['bins'])
  write_dx(gist_data['Esw_norm'],       gist_data['origin'], gist_data['dim'], '%s_Esw_norm'       %output, gist_data['bins'])
  write_dx(gist_data['Eww_dens'],       gist_data['origin'], gist_data['dim'], '%s_Eww_dens'       %output, gist_data['bins'])
  write_dx(gist_data['Eww_norm_unref'], gist_data['origin'], gist_data['dim'], '%s_Eww_norm_unref' %output, gist_data['bins'])
  write_dx(Eww_norm,                    gist_data['origin'], gist_data['dim'], '%s_Eww_norm_ref'   %output, gist_data['bins'])
  write_dx(gist_data['Dipole_x_dens'],  gist_data['origin'], gist_data['dim'], '%s_Dipole_x_dens'  %output, gist_data['bins'])
  write_dx(gist_data['Dipole_y_dens'],  gist_data['origin'], gist_data['dim'], '%s_Dipole_y_dens'  %output, gist_data['bins'])
  write_dx(gist_data['Dipole_z_dens'],  gist_data['origin'], gist_data['dim'], '%s_Dipole_z_dens'  %output, gist_data['bins'])
  write_dx(gist_data['Dipole_dens'],    gist_data['origin'], gist_data['dim'], '%s_Dipole_dens'    %output, gist_data['bins'])
  write_dx(gist_data['neighbor_dens'],  gist_data['origin'], gist_data['dim'], '%s_neighbor_dens'  %output, gist_data['bins'])
  write_dx(gist_data['neighbor_norm'],  gist_data['origin'], gist_data['dim'], '%s_neighbor_norm'  %output, gist_data['bins'])
  write_dx(gist_data['order_norm'],     gist_data['origin'], gist_data['dim'], '%s_order_norm'     %output, gist_data['bins'])


def main():

  ############################# M A I N   P A R T  --  Command line arguments ###############################

  time1 = time.time()

  parser.add_argument("--input",               "-i",        help="Input data path, format must be dx or xplor specified with --format flag.", type=argparse.FileType('r'))
  parser.add_argument("--input_format",        "-fi",       help="Input data format, must be 'dx', 'xplor' or 'pdb'", choices=['dx','xplor', 'pdb', 'gist_out'])
  parser.add_argument("--output",              "-o",        help="Output file name, format is always dx or pdb.", type=str, default="output")
  parser.add_argument("--job",                 "-j",        help="Specify the job.", choices=['smooth_gauss','gist_sphere_expander','gist_cavity_expander','population_analysis_michael','sd_smooth','map_gist','gist_cavity','is_in_box','gist_grid_analysis', \
                                                                                              'write_dx_boltzmann_invert', 'gist_water_analysis', 'find_random','nearest','placevent','find_max_gradient_preselection','diff', 'find_max', \
                                                                                              'write_dx', 'write_pdb', 'real_cut', 'dens_v_t', 'write_dx_sd', 'find_max_dc', 'pdb2table', 'map2table', 'site_site', 'find_max_grad', 'smooth_gauss', \
                                                                                              'calc_rdd'], required=True)
  parser.add_argument("--ref_map",             "-rmap",     help="Reference map for '--job diff'. Calculation via 'Input data' - 'Reference map' = 'output'.", type=argparse.FileType('r'))
  parser.add_argument("--ref_model",           "-rmod",     help="Reference model for '--job calc_rdd'. Must be pdb-format.", type=argparse.FileType('r'))
  parser.add_argument("--ref_model_selection", "-rmod_sel", help="Selection for reference model (--ref_model) in charmm-syntax. Use it for '--job calc_rdd'.", type=str)
  parser.add_argument("--ref_map_format",      "-frmap",    help="Reference map file format. Must be format dx or xplor.")
  parser.add_argument("--threshold",           "-th",       help="Threshold for lattice points that should be considered in the different methods for maximum calculation. No default!", type=float)
  parser.add_argument("--trajectory",          "-tj",       help="Input amber trajectory in netcdf format", type=argparse.FileType('r'))
  parser.add_argument("--parameter",           "-pa",       help="Input amber parameter and topology file.", type=argparse.FileType('r'))
  parser.add_argument("--hsite",               "-hi",       help="Input hydration sites as pdb format.", type=argparse.FileType('r'))
  parser.add_argument("--hradius",             "-hr",       help="Radius of hydration sites. Default is 1.0 Angstrom.", type=float, default=1.0)
  parser.add_argument("--numproc",             "-np",       help="Number of processors used. Default = 1. Just for dens_v_t method.", type=int, default=2)
  parser.add_argument("--start",               "-start",    help="Start frame for MD Analysis. Default = 0. Just for dens_v_t method.", type=int, default=0)
  parser.add_argument("--stop",                "-stop",     help="Stop frame for MD Analysis. Default = -1 (last frame). Just for dens_v_t method.", type=int, default=-1)
  parser.add_argument("--skip",                "-skip",     help="Frames skipped during MD Analysis. Default = 1 (read all frames ). Just for dens_v_t method.", type=int, default=1)
  parser.add_argument("--grid_factor",         "-gf",       help="Factor that multiplies each grid value during write_dx process. Default=1.0.", type=float, default=1.0)
  parser.add_argument("--random_sites",        "-rs",       help="How many random sites are generated. No default", type=int)
  parser.add_argument("--eww_matrix",          "-eww",      help="Water-water interaction matrix as output from cpptraj's gist approach.", type=argparse.FileType('r'))
  parser.add_argument("--reference_energy",    "-re",       help="Water-water interaction energy in bulk-phase for appropriate water model. Default -11.036 kcal/mol for TIP4PEW.", type=float, default=-11.036)
  parser.add_argument("--reference_density",   "-rd",       help="Water number-density in bulk-phase for appropriate water model. Default 0.0332 Mol/A^-3 for TIP4PEW. If not given, it is assumed that it equals calculation density.", type=float)
  parser.add_argument("--calculation_density", "-cd",       help="Water number-density in bulk-phase That was used for calculation of g(r). Default 0.0332 Mol/A^-3 for TIP4PEW.", type=float, default=0.0332)
  parser.add_argument("--vdw_add",             "-vdw_add",  help="Additive constant that increases the atoms vdw radii before msms calculation. Default=0.0", type=float, default=0.0)
  parser.add_argument("--vdw_add_max",         "-vdw_max",  help="Maximum additive constant to vdw-radius for gist_cavity_expander. Default = 4.0", type=float, default=4.0)
  parser.add_argument("--vdw_add_dr",          "-vdw_dr",   help="Incremental for radius expanision in gist_cavity_expander. Default = 0.25", type=float, default=0.25)
  parser.add_argument("--happy_only",          "-happy",    help="Use happy waters only. Default = False", type=bool, default=False)
  parser.add_argument("--unhappy_only",        "-unhappy",  help="Use unhappy waters only. Default = False", type=bool, default=False)
  args = parser.parse_args()

  if args.reference_density == None:
    args.reference_density = args.calculation_density

  ############################# M A I N   P A R T  --  Command line control ###############################

  if args.job == 'gist_sphere_expander':

    gist_data = read_map (args.input.name, args.input_format)
    if not args.eww_matrix == None:
      Eww_matrix = read_map (args.eww_matrix.name, 'Eww_ij', False, gist_data, args.calculation_density)
    else:
      Eww_matrix = {}
    water_positions = read_pdb(args.hsite.name)
    gist_sphere_expander(water_positions, gist_data, Eww_matrix, args.reference_energy, args.calculation_density, args.reference_density, args.output, args.vdw_add_max, args.vdw_add_dr)



  if args.job == 'gist_cavity_expander':
  
    gist_data = read_map (args.input.name, args.input_format)
  #  if not args.eww_matrix == None:
  #    Eww_matrix = read_map (args.eww_matrix.name, 'Eww_ij', False, gist_data, args.calculation_density)
  #  else:
  #    Eww_matrix = {}
    if args.eww_matrix == None:
      eww_matrix_path = None
    else:
      eww_matrix_path = args.eww_matrix.name
  
    structure_path = args.ref_model.name
    gist_cavity_expander(structure_path, gist_data, eww_matrix_path, args.reference_energy, args.calculation_density, args.reference_density, args.output, args.vdw_add_max, args.vdw_add_dr, args.happy_only, args.unhappy_only)
  
  
  if args.job == 'population_analysis_michael':
  
    map_1 = read_map (args.input.name, args.input_format)
    water_positions = read_pdb(args.hsite.name)
    population_analysis_michael (water_positions, map_1, args.hradius, args.output)
  
  if args.job == 'gist_grid_analysis':
  
    gist_data = read_map (args.input.name, args.input_format)
    #if not args.eww_matrix == None:
      #Eww_matrix = read_map (args.eww_matrix.name, 'Eww_ij', False, gist_data, args.calculation_density)
    #else:
      #Eww_matrix = {}
    gist_grid_analysis(gist_data, args.eww_matrix.name, args.calculation_density, args.reference_density, args.reference_energy)
  
  
  if args.job == 'gist_water_analysis':
  
    gist_data = read_map (args.input.name, args.input_format)
    if not args.eww_matrix == None:
      Eww_matrix = read_map (args.eww_matrix.name, 'Eww_ij', False, gist_data, args.calculation_density)
    else:
      Eww_matrix = {}
    water_positions = read_pdb(args.hsite.name)
    gist_water_analysis (water_positions, gist_data, args.hradius, Eww_matrix, args.reference_energy, args.calculation_density, args.reference_density, args.output)
  
  if args.job == 'gist_cavity':
  
    gist_data = read_map (args.input.name, args.input_format)
  #  if not args.eww_matrix == None:
  #    Eww_matrix = read_map (args.eww_matrix.name, 'Eww_ij', False, gist_data, args.calculation_density)
  #  else:
  #    Eww_matrix = {}
    if args.eww_matrix == None:
      eww_matrix_path = None
    else:
      eww_matrix_path = args.eww_matrix.name
  
    structure_path = args.ref_model.name
    gist_cavity(structure_path, gist_data, eww_matrix_path, args.reference_energy, args.calculation_density, args.reference_density, args.output, args.vdw_add, args.happy_only, args.unhappy_only)
  
  
  if args.job == 'map_gist':
  
    gist_data = read_map (args.input.name, args.input_format)
    if not args.eww_matrix == None:
      Eww_matrix = read_map (args.eww_matrix.name, 'Eww_ij', False, gist_data, args.calculation_density)
    else:
      Eww_matrix = {}
    map_gist(gist_data, Eww_matrix, args.reference_energy, args.calculation_density, args.reference_density, args.output)
  
  
  if args.job == 'is_in_box':
  
    pdb_crd = read_pdb(args.ref_model.name)
    box = read_map(args.input.name, args.input_format, True)
  
    who_is_in_the_box(pdb_crd, box)
  
  
  if args.job == 'write_dx_boltzmann_invert':
  
    map1 = read_map (args.input.name, args.input_format)
    boltzmann_inverted_map = boltzmann_inverter(map1['map'], map1['bins'])
    write_dx(boltzmann_inverted_map, map1['origin'], map1['dim'], args.output, map1['bins'])
  
  
  if args.job == 'find_random':
  
    map1 = read_map (args.input.name, args.input_format)
    random_map = find_random(map1, args.random_sites)
    write_pdb(random_map, args.output)
    
  
  if args.job == 'nearest':
  
    sites = read_pdb (args.input.name)
    ref_sites = read_pdb (args.ref_model.name)
    map_data = read_map (args.ref_map.name, args.ref_map_format, True)
    data = nearest_neighbour_dist(sites, ref_sites, map_data, args.threshold)
    write_table(data, args.output)
  
    o = open('%s_crd_len.dat' %args.output, "w")
    o.write("Nr. hydration sites\n%d" %len(sites))
    o.close()
  
  
  if args.job == 'placevent':
  
    map1 = read_map (args.input.name, args.input_format)
    maximum = placevent(map1, args.output)
  
  
  if args.job == 'find_max_gradient_preselection':
  
    map1 = read_map (args.input.name, args.input_format)
    maximum = find_max_gradient_preselection(map1, args.hradius)
    write_pdb(maximum, args.output)
  
  
  if args.job == 'real_cut':
  
    map1 = read_map (args.input.name, args.input_format)
    ref_map = read_map (args.ref_map.name, args.ref_map_format, True)
    cutted_map = real_cut (map1, ref_map, args.output)
    write_dx(cutted_map, ref_map['origin'], ref_map['dim'], args.output, ref_map['bins'])
  
  
  if args.job =='calc_rdd':
  
    map1 = read_map (args.input.name, args.input_format)
    calc_radial_density_distribution (map1, args.ref_model.name, args.output,  args.ref_model_selection, args.threshold)
  
  
  if args.job == 'sd_smooth':
  
    map1 = read_map (args.input.name, args.input_format)
    smooth_map = smooth_gauss(map1['map'], map1['bins'])
    sd_map = calc_sd(smooth_map)
    write_dx(sd_map, map1['origin'], map1['dim'], args.output, map1['bins'], factor=args.grid_factor)
  
  if args.job == 'smooth_gauss':
  
    map1 = read_map (args.input.name, args.input_format)
    smooth_map = smooth_gauss(map1['map'], map1['bins'])
    write_dx(smooth_map, map1['origin'], map1['dim'], args.output, map1['bins'], factor=args.grid_factor)
  
  if args.job == 'find_max_grad':
  
    map1 = read_map (args.input.name, args.input_format)
    max_map = find_max_grad(map1)
    write_pdb(max_map, args.output)
  
  
  if args.job == 'site_site':
  
    map1 = {}
  
    if (args.input_format == 'pdb'):
  
      map1 = read_pdb (args.input.name)
  
    if (args.input_format == 'xplor' or args.input_format == 'dx'):
  
      map1 = read_map (args.input.name, args.input_format)
  
    output = site_site (map1, args.output)
  
  
  if args.job == 'pdb2table':
  
    map1 = read_pdb (args.input.name)
    write_table(map1, args.output)
  
  
  if args.job == 'map2table':
  
    map1 = read_map (args.input.name, args.input_format)
    write_map_table(map1, args.output)
  
  
  if args.job == 'write_dx_sd':
  
    map1 = read_map (args.input.name, args.input_format)
    sd_map = calc_sd(map1['map'])
    write_dx(sd_map, map1['origin'], map1['dim'], args.output, map1['bins'], factor=args.grid_factor)
  
  
  if args.job == 'diff':
  
    map1 = read_map (args.input.name, args.input_format)
    map2 = read_map (args.ref_map.name, args.ref_map_format)
  
    if map1['bins'][0] != map2['bins'][0] or map1['bins'][1] != map2['bins'][1] or map1['bins'][2] != map2['bins'][2] or map1['origin'][0] != map2['origin'][0] or map1['origin'][1] != map2['origin'][1] or map1['origin'][2] != map2['origin'][2] or map1['n'][0] != map2['n'][0] or map1['n'][1] != map2['n'][1] or map1['n'][2] != map2['n'][2]:
  
      print '%s and %s have to match exactly!' %(args.input, args.ref_map)
      exit()
  
    diff = map_diff (map1['map'], map2['map'],map1['bins'], )
    write_dx(diff, map1['origin'], map1['dim'], args.output, map1['bins'], factor=args.grid_factor)
  
  if args.job == 'dens_v_t':
  
    if args.output == 'output':
      print "Specify output-name."
      exit()
    
    density_vs_time (args.trajectory.name, args.parameter.name, args.hsite.name, args.hradius, args.output, args.numproc, args.start, args.stop, args.skip)
  
  
  if args.job == 'find_max':
  
    map1 = read_map (args.input.name, args.input_format)
    maximum = find_maximum(map1, args.threshold, args.hradius)
    write_pdb(maximum, args.output)
  
  
  if args.job == 'find_max_dc':
  
    map1 = read_map (args.input.name, args.input_format)
    real_map = frac2real (map1['map'], map1['origin'], map1['dim'], map1['n'], map1['alpha'], map1['beta'], map1['gamma'])
    maximum = find_max_dc (real_map, args.hradius, args.threshold)
    write_pdb(maximum, args.output)
  
  
  if args.job == 'write_pdb':
  
    map1 = read_map (args.input.name, args.input_format)
    write_pdb(frac2real(map1['map'], map1['origin'], map1['dim'], map1['n'], map1['alpha'], map1['beta'], map1['gamma']), args.output)
  
  
  if args.job == 'write_dx':
  
    map1 = read_map (args.input.name, args.input_format)
    write_dx(map1['map'], map1['origin'], map1['dim'], args.output, map1['bins'], factor=args.grid_factor)
  
  time2 = time.time()
  print '\n Calculation took %.2f seconds.' %(time2-time1)
  print '\n Calculation terminated at %s.' %time.ctime()

if __name__ == "__main__":

  main()