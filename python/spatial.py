### Written by T. Wulsdorf @ AG Klebe Marburg University
### Please contact via tobias.wulsdorf@uni-marburg.de

import numpy as np

def rotate_check(matrix):

	if not (0.99 < np.linalg.det(matrix) < 1.01):

		raise Warning("Warning: Determinant of rotation matrix is %s. Should be close to +1.0." %np.linalg.det(matrix))


def do_rotation (crds, origin, rot_mat): 

	return (crds - origin).dot(rot_mat) + origin

#this is an alias
def rotate (crds, origin, rot_mat):

	return do_rotation(crds, origin, rot_mat)

def make_grid(arrays, out=None):
	"""
	!!! Adapted from:
	!!! http://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays

	Generate a cartesian product of input arrays.

	Parameters
	----------
	arrays : list of array-like
		1-D arrays to form the cartesian product of.
	out : ndarray
		Array to place the cartesian product in.

	Returns
	-------
	out : ndarray
		2-D array of shape (M, len(arrays)) containing cartesian products
		formed of input arrays.

	Examples
	--------
	>>> make_grid(([1, 2, 3], [4, 5], [6, 7]))
	array([[1, 4, 6],
			[1, 4, 7],
			[1, 5, 6],
			[1, 5, 7],
			[2, 4, 6],
			[2, 4, 7],
			[2, 5, 6],
			[2, 5, 7],
			[3, 4, 6],
			[3, 4, 7],
			[3, 5, 6],
			[3, 5, 7]])

	"""

	arrays = [np.asarray(x) for x in arrays]

	dtype  = arrays[0].dtype

	n = np.prod([x.size for x in arrays])

	if out is None:

		out = np.zeros([n, len(arrays)], dtype=dtype)

	m = n / arrays[0].size

	out[:,0] = np.repeat(arrays[0], m)

	if arrays[1:]:

		make_grid(arrays[1:], out=out[0:m,1:])

		for j in xrange(1, arrays[0].size):

			out[j*m:(j+1)*m,1:] = out[0:m,1:]

	return out

