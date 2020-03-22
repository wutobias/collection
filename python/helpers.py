### Written by T. Wulsdorf @ AG Klebe Marburg University
### Please contact via tobias.wulsdorf@uni-marburg.de

import numpy as np

def are_you_numpy(a):

	"""
	Returns True if a is an instance of numpy.
	False otherwise.
	"""

	return type(a).__module__ == np.__name__