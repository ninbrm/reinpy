import h5py
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt
from grid_manager import Grid

FPATH = '../data/affinity.h5' # Replace with the correct path

class HDF5Reader:

	def __init__(self, fpath=FPATH):
		self.fpath = fpath

		f = h5py.File(fpath)
		# self.ncols = f['NCOLS'][0]
		# self.nrows = f['NROWS'][0]
		self.nrows = f['NCOLS'][0] # rows and cols seems to be mixed up
		self.ncols = f['NROWS'][0]
		self.xllcorner = f['XLLCORNER'][0] #the longitude coordinates of the lower left corner
		self.yllcorner = f['YLLCORNER'][0] # the latitude coordinates of the lower left corner
		self.xcellsize = f['XCELLSIZE'][0] # width of the pixels
		self.ycellsize = f['YCELLSIZE'][0] # height of the pixels (when pixels are square this will be the  same value as the XCELLSIZE, i.e. most cases)
		self.epsg = f['EPSG'][0]  # is the code to represent the projection
		self.I = f['AFFINITY_START'][:].tolist()
		self.J = f['AFFINITY_END'][:].tolist()
		self.V = f['AFFINITY_VALUES'][:].tolist()

		# Below are hacks that deal with the missing node in the data:
		self.I = [i-1 for i in self.I]
		self.I.append(max(self.I) + 1)
		self.J = [j-1 for j in self.J]
		self.J.append(max(self.J) + 1)
		self.V.append(0)

		# self.V = [v-0.9 for v in self.V]
		Q = f['QUALITY_VALUES'][:].tolist()
		self.qualities = np.array(Q)
		self.nodata_value = f['NODATA_VALUE'][0]  # what represents the missing data in the QUALITY vector
		self.npix = max(max(self.I),max(self.J)) # number of pixels

	def truncate_affinities(self, min_affinity=-4, max_affinity=4):
		'''Truncate affinity values between min_affinity and max_affinity
		'''
		self.V = np.fmin(max_affinity, self.V)
		self.V = np.fmax(min_affinity, self.V)

	def standardize_affinities(self):
		'''Remove mean and normalize std to 1:
		'''
		meanV = np.mean(self.V)
		stdV = np.std(self.V)
		self.V = (self.V-meanV)/stdV

	def make_grid(self):
		# Exponential of the affinities:
		A = sparse.coo_matrix((np.exp(self.V), (self.I,self.J))).tocsr() # build a graph with the exponential of the affinities as edge weight

		# Uniform affinities:
		# A = sparse.coo_matrix((np.ones(len(self.V)), (self.I,self.J))).tocsr() # build a graph with the exponential of the affinities as edge weight

		# Non-exponential affinities:
		# V_min = min(self.V)
		# self.V = [v - V_min + 1 for v in self.V]
		# A = sparse.coo_matrix((np.array(self.V), (self.I,self.J))).tocsr() # build a graph with the exponential of the affinities as edge weight

		map_size = (self.nrows*self.ycellsize, self.ncols*self.xcellsize)
		return Grid(map_size, (self.nrows, self.ncols), graph=A, qualities=self.qualities)