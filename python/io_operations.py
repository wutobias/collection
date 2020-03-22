### Written by T. Wulsdorf @ AG Klebe Marburg University
### Please contact via tobias.wulsdorf@uni-marburg.de

import numpy as np
from spatial import rotate_check, rotate, do_rotation, make_grid
from helpers import are_you_numpy
from string import ascii_uppercase
import gzip

def load_dx(path):

    if path.endswith(".gz"):
        map_file_tmp = gzip.open(path,"r")
    else:
        map_file_tmp = open(path,"r")
    map_file     = map_file_tmp.readlines()

    origin    = np.zeros(3, dtype=float)
    bins      = np.zeros(3, dtype=int)
    frac2real = np.eye(3,3)

    start_row = -1

    for i, item in enumerate(map_file):

        if len(item.rstrip().split()) > 1 and item.rstrip().split()[0] == 'object' and item.rstrip().split()[1] == '1':

            start_row = i

    origin[0] = float(map_file[start_row + 1].rstrip().split()[1])
    origin[1] = float(map_file[start_row + 1].rstrip().split()[2])
    origin[2] = float(map_file[start_row + 1].rstrip().split()[3])

    bins[0]   = int(map_file[start_row].rstrip().split()[5])
    bins[1]   = int(map_file[start_row].rstrip().split()[6])
    bins[2]   = int(map_file[start_row].rstrip().split()[7])

    frac2real[0][0] = float(map_file[start_row+2].rstrip().split()[1])
    frac2real[0][1] = float(map_file[start_row+2].rstrip().split()[2])
    frac2real[0][2] = float(map_file[start_row+2].rstrip().split()[3])
    frac2real[1][0] = float(map_file[start_row+3].rstrip().split()[1])
    frac2real[1][1] = float(map_file[start_row+3].rstrip().split()[2])
    frac2real[1][2] = float(map_file[start_row+3].rstrip().split()[3])
    frac2real[2][0] = float(map_file[start_row+4].rstrip().split()[1])
    frac2real[2][1] = float(map_file[start_row+4].rstrip().split()[2])
    frac2real[2][2] = float(map_file[start_row+4].rstrip().split()[3])

    f = field(Bins=bins, Frac2Real=frac2real, Origin=origin)
    v = np.zeros(bins, dtype=float)

    rows_of_data = (bins[0] * bins[1] * bins[2]) / 3

    if (bins[0] * bins[1] * bins[2]) % 3 != 0:
        rows_of_data = rows_of_data + 1

    v_i = 0
    for i in range(start_row+7, start_row+7+rows_of_data):
        for value in map_file[i].rstrip().split():

            x      =   v_i / (bins[1]*bins[2])
            y      = ( v_i % (bins[1]*bins[2])) / bins[1]
            z      = ( v_i % (bins[1]*bins[2])) % bins[2]

            v[x,y,z] = value
            v_i += 1

    return f,v


class field(object):

    __doc__="""
    This is a class for operation on generic
    scalar fields. Like GIST objects, they
    must be described in cartesian as well as 
    in fractional space. That means that we need 
    origin, frac2real/real2frac matrix(vector)
    and vector.

    Tobias Wulsdorf, AG Klebe, 08/2016
    """

    def __init__(self, Bins, Frac2Real=None, Delta=None, Origin=None, Center=None):


        if type(Frac2Real) == type(None) and type(Delta) == type(None):

            raise ValueError("Frac2Real or Delta must be given!")

        if type(Frac2Real) != type(None) and type(Delta) != type(None):

            raise ValueError("Either Frac2Real or Delta must be given!")

        if type(Frac2Real) == type(None):

            self.delta     = Delta
            self.frac2real = np.eye(3,3) * self.delta

        else:

            self.frac2real = Frac2Real
            self.delta     = np.linalg.norm(self.frac2real, axis=0)

        self.real2frac = np.linalg.inv(self.frac2real)
        self.bins      = Bins

        self.rotation_matrix    = np.eye(3,3)
        self.translation_vector = np.zeros(3)


        if type(Origin) == type(None) and type(Center) == type(None):

            raise ValueError("Origin or Center must be given!")

        if type(Origin) != type(None) and type(Center) != type(None):

            raise ValueError("Either origin or center must be given!")

        if type(Center) == type(None):

            self.origin    = Origin
            self.center    = self.get_real(self.bins/2)

        else:

            self.center    = Center
            #First we need an auxiliary origin at (0,0,0)
            self.origin    = np.zeros(3)
            #Second translate origin according center displacement
            self.origin    = self.center - self.get_real(self.bins/2)


    def translate(self, vector=np.zeros(3)):

        __doc__="""
Translatation vector of unit cell origin
"""

        self.translation_vector += vector


    def rotate(self, matrix=np.eye(3,3)):

        __doc__=""" 
Rotate the unit cell vectors. 
Does not affect the gridData objects.
"""

        rotate_check(matrix)
        self.rotation_matrix = matrix.dot(self.rotation_matrix)


    def translate_global(self, vector=np.zeros(3)):

        __doc__="""
Translate global coordinate system
along vector.
"""

        self.origin += vector


    def rotate_global(self, reference_point=np.zeros(3), matrix=np.eye(3,3)):

        __doc__="""
Rotate global coordinate system around
reference point.
"""

        rotate_check(matrix)
        self.origin = do_rotation(self.origin, reference_point, matrix)

        self.rotate(matrix)

        self.translation_vector = do_rotation(self.translation_vector, np.zeros(3), matrix)


    def get_nice_frac2real(self):

        return self.rotation_matrix.dot(self.frac2real)


    def get_nice_real2frac(self):

        return np.linalg.inv(self.get_nice_frac2real())


    def get_voxel_volume(self):

        __doc__="""
Returns the volume per grid voxel.
"""

        return np.absolute(np.cross(self.frac2real[:,0], self.frac2real[:,1]).dot(self.frac2real[:,2]))


    def get_frac(self, real_array):

        #Convert to initial real space by inverse translation and rotation
        initial_reals = do_rotation(real_array, self.origin + self.translation_vector, np.linalg.inv(self.rotation_matrix))

        #Remove origin
        initial_reals -= (self.origin + self.translation_vector)

        #Convert to initial fractional space
        return initial_reals.dot(self.real2frac)


    def get_real(self, frac_array):

        #Convert to real space
        reals = np.array(frac_array).dot(self.frac2real)

        #Perform rotation translation
        return do_rotation(reals, np.zeros(3), self.rotation_matrix) + self.origin + self.translation_vector


    def get_centers(self):

        return self.get_real(make_grid((np.arange(self.bins[0]),\
                                        np.arange(self.bins[1]),\
                                        np.arange(self.bins[2]))))

    def get_centers_real(self):

        return self.get_centers()


    def get_centers_frac(self):

        return make_grid((np.arange(self.bins[0]),\
                          np.arange(self.bins[1]),\
                          np.arange(self.bins[2])))


class write_files(object):

    def __init__(self, Delta=None, Frac2Real=None, Bins=None, Origin=None, Value=None, XYZ=None, X=None, Y=None, Z=None, Format='PDB', Filename=None, Remarks=None, Nan_fill=-1.0):

        """
        This class can write different file types.
        currently only dx and pdb are supported.
        """

        self._delta     = Delta
        self._frac2real = Frac2Real
        self._bins      = Bins
        self._origin    = Origin
        self._value     = Value
        self._x         = X
        self._y         = Y
        self._z         = Z
        self._format    = Format
        self._filename  = Filename
        self._xyz       = XYZ
        self._remarks   = Remarks
        self._nan_fill  = Nan_fill

        if type(self._filename) != str:

            self._filename  = 'output.'
            self._filename  += self._format

        self._writers = {
                        'PDB'  : self._write_PDB,
                        'DX'   : self._write_DX,
                        'GIST' : self._write_GIST
                        }

        data = self._writers[self._format]()

        o = open(self._filename, "w")
        o.write(data)
        o.close()

    def _merge_x_y_z(self):

        return np.stack( ( self._x, self._y, self._z ), axis=1 )


    def _write_PDB(self):

        """
        Write a PDB file.
        This is intended for debugging. It writes all atoms
        as HETATM of element X with resname MAP.
        """

        if are_you_numpy(self._xyz):

            if self._xyz.shape[-1] != 3:

                raise TypeError(
                    "XYZ array has wrong shape.")

        else:

            if not ( are_you_numpy(self._x) or are_you_numpy(self._y) or are_you_numpy(self._z) ):

                raise TypeError(
                    "If XYZ is not given, x,y and z coordinates must be given in separate arrays.")

            else:

                self._xyz = self._merge_x_y_z()

        if type(self._value) == type(None):

            self._value = np.zeros( len(self._xyz), dtype=float )

        data  = 'REMARK   6\n'
        data += 'REMARK   6 File written by write_files.py\n'
        if self._remarks != None:

            for line_idx in range(0, len(self._remarks), 86):
                data += 'REMARK   6 %s\n' %self._remarks[line_idx:line_idx+86]

        for xyz_i, xyz in enumerate(self._xyz):

            #iterate over uppercase letters
            chain_id    = ascii_uppercase[( len(str(xyz_i+1)) / 5 )]

            atom_counts = xyz_i - ( len(str(xyz_i+1)) / 6 ) * 100000
            resi_counts = xyz_i - ( len(str(xyz_i+1)) / 5 ) * 10000
            data += \
            '%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          \n' \
            %('HETATM',atom_counts+1,'X','', 'MAP', chain_id, resi_counts+1, '', xyz[0], xyz[1], xyz[2], 0.00, float( self._value[xyz_i] ) )

        data += 'END\n'

        return data


    def _write_DX(self):

        """
        Writes DX files according to openDX.
        """

        if not ( are_you_numpy(self._origin) or are_you_numpy(self._bins) ):

            raise TypeError(
            "Origin and bins must be given.")

        #This means not (a XOR b) or not (a or b)
        if are_you_numpy(self._delta) == are_you_numpy(self._frac2real) :

            raise TypeError(
            "Either delta or frac2real must be given.")

        if are_you_numpy(self._delta):

            self._frac2real = np.zeros((3,3), dtype=float)

            np.fill_diagonal(self._frac2real, self._delta)

        data = '''object 1 class gridpositions counts %d %d %d
origin %8.4f %8.4f %8.4f
delta %8.4f %8.4f %8.4f
delta %8.4f %8.4f %8.4f
delta %8.4f %8.4f %8.4f
object 2 class gridconnections counts %d %d %d
object 3 class array type float rank 0 items %d data follows
''' %(self._bins[0], self._bins[1], self._bins[2],\
      self._origin[0], self._origin[1], self._origin[2],\
      self._frac2real[0][0], self._frac2real[0][1], self._frac2real[0][2],\
      self._frac2real[1][0], self._frac2real[1][1], self._frac2real[1][2],\
      self._frac2real[2][0], self._frac2real[2][1], self._frac2real[2][2],\
      self._bins[0],   self._bins[1],   self._bins[2],\
      self._bins[2] * self._bins[1] * self._bins[0])

        i = 0
        for x_i in range(0, self._bins[0]):

            for y_i in range(0, self._bins[1]):

                for z_i in range(0, self._bins[2]):

                    ### writing an integer instead of float
                    ### saves us some disk space
                    if np.isnan(self._value[x_i][y_i][z_i]):

                        data += str(self._nan_fill) + " "

                    else:

                        if self._value[x_i][y_i][z_i] == 0.0:

                            data += "0 " 

                        else:

                            data += str(self._value[x_i][y_i][z_i]) + ' '
                
                    i += 1

                    if i == 3:

                        data += '\n'
                        i = 0
        return data

    def _write_GIST(self):

        pass

class PDB(object):

  def __init__(self, Path):

    self.path = Path

    self.crd  = list()
    self.B    = list()

    with open(self.path, "r") as PDB_file:

      for i, line in enumerate(PDB_file):

        if not (line[0:6].rstrip() == 'ATOM' or line[0:6].rstrip() == 'HETATM'):

          continue

        if i <= 9999:

          #Coordinates
          self.crd.append(list())
          self.crd[-1].append(float(line.rstrip()[30:38]))
          self.crd[-1].append(float(line.rstrip()[38:46]))
          self.crd[-1].append(float(line.rstrip()[46:54]))

          #B-Factors
          self.B.append(line.rstrip()[54:59])

        if 9999 < i <= 99999:

          #Coordinates
          self.crd.append(list())
          self.crd[-1].append(float(line.rstrip()[31:39]))
          self.crd[-1].append(float(line.rstrip()[39:47]))
          self.crd[-1].append(float(line.rstrip()[47:55]))

          #B-Factors
          self.B.append(line.rstrip()[55:60])

        if i > 99999:

          #Coordinates
          self.crd.append(list())
          self.crd[-1].append(float(line.rstrip()[33:41]))
          self.crd[-1].append(float(line.rstrip()[41:49]))
          self.crd[-1].append(float(line.rstrip()[49:57]))

          #B-Factors
          self.B.append(line.rstrip()[57:62])

    self.crd  = np.array(self.crd)
    self.B    = np.array(self.B)