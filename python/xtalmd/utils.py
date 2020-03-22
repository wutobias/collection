import numpy as np
from xtalmd.constants import DEG2RAD, RAD2DEG

def cellbasis(angles, edges):
	'''
	For the unit cell with given angles and edge lengths calculate the basis
	transformation (vectors) as a 4x4 np.array
	'''
	rad = [DEG2RAD*i for i in angles]
	basis = np.identity(4)
	basis[0][1] = np.cos(rad[2])
	basis[1][1] = np.sin(rad[2])
	basis[0][2] = np.cos(rad[1])
	basis[1][2] = (np.cos(rad[0]) - basis[0][1]*basis[0][2])/basis[1][1]
	basis[2][2] = np.sqrt(1 - basis[0][2]**2 - basis[1][2]**2)
	edges.append(1.0)
	return basis * edges # np.array multiplication!
