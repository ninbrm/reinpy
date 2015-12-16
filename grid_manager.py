from __future__ import division
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import sparse
from bisect import bisect_left

class Grid:
	'''
	Manage landscapes as grids, represented by a sparse matrix of affinities.
	Allows to do plots and to do conversions between the different coordinates of nodes (id in the graph, physical coordinates, grid coordinates)

	columns of the grid are parallels to the y axis of the landscape, and rows to the x axis.
	rows and columns are numbered from the top left corner.
	The order of nodes in the graph correspond to reading the grid from left to right, top to bottom.
	i.e. if there are n columns, the pixel on the ith rows, jth column correspond to the node i*n+j
	All numbering start with 0.

	Physical coordinates take the form (value on the y axis, value on the x axis)
	Grid coordinates take the form (row id, column id)

	Instance variables:
		- A : affinity matrix in the csr format. Great for matrix-vector multiplications
		- A_dok : affinity matrix in the dok format. Access to edge weight in O(1)
		- map_size : tuple (size along the y axis, size along the x axis). This is the physical size (in meters or kilometers) of the landscape.
		- shape : tuple (number of rows, number of columns). This is the size of the grid representing the landscape
		- pixel_size : tuple (pixel size along the y axis, pixel size along the x axis)
		- x_ticks : list delimitating the pixels along the x axis (in the physical coordinates). The ith column is delimitated by x_ticks[i] and x_ticks[i+1]. 
		if there are n colums, len(x_ticks) will be n+1.
		- y_ticks : list delimitating the pixels along the y axis (in the physical coordinates). The ith row is delimitated by y_ticks[i] and y_ticks[i+1]. 
		if there are n rows, len(y_ticks) will be n+1.
	'''

	def __init__(self, map_size, shape, graph=None, qualities=None):
		'''
		Inputs: 
			- graph : sparse matrix (if None, the graph can be set later using set_affinities)
			- map_size : tuple (size along the y axis, size along the x axis). This is the physical size (in meters or kilometers) of the landscape.
			- shape : tuple (number of rows, number of columns). This is the size of the grid representing the landscape
		'''
		self.pixel_size = (map_size[0]/shape[0], map_size[1]/shape[1])
		self.shape = shape
		self.map_size = map_size
		self.n_rows, self.n_cols = shape
		self.N = shape[0] * shape[1]
		self.x_ticks = np.linspace(0, map_size[1], num=(shape[1]+1))
		self.y_ticks = np.linspace(0, map_size[0], num=(shape[0]+1))
		
		if graph is not None:
			self.set_affinities(graph)

		if qualities is not None:
			self.qualities = qualities
		else:
			self.qualities = np.ones(self.N)


	def set_affinities(self, graph):
		'''
		Set the affinity matrix A and A_dok.
		graph should be a sparse matrix.
		'''
		if graph.shape[0] != self.N:
			raise ValueError('The shape of the grid does not match the number of nodes of the graph')
		self.A = graph.tocsr()
		self.A_dok = graph.todok()
	
	def node_id_to_grid_coordinates(self, node_id):
		'''
		Return the grid coordinates (row id, column id) corresponding to a given node
		'''
		self.check_node_id(node_id)
		j = int(node_id % self.n_cols)
		i = (node_id - j) / self.n_cols
		return (i, j)


	def node_id_to_coordinates(self, node_id):
		'''
		Return the physical coordinates of a node (center of the corresponding pixel)
		'''
		return self.grid_coordinates_to_coordinates(self.node_id_to_grid_coordinates(node_id))
	
	def coordinates_to_grid_coordinates(self, coo):
		'''
		Return the grid coordinates (row id, column id) of the pixel containing the 
		physical coordinates "coo" (value along the y axis, value along the x axis)
		'''
		self.check_coordinates(coo)
	
		i = max(0, bisect_left(self.y_ticks, coo[0]) - 1)
		j = max(0, bisect_left(self.x_ticks, coo[1]) - 1)
		return (i,j)

	def coordinates_to_node_id(self, coo):
		'''
		Return the id of the node containing the physical coordinates "coo"
		the coordinates are interpreted as (value along the y axis, value along the x axis)
		'''
		return self.grid_coordinates_to_node_id(self.coordinates_to_grid_coordinates(coo))

	def grid_coordinates_to_node_id(self, grid_coo):
		'''
		Return the id of the node corresponding to the grid coordinates "grid_coo" (row id, column id)
		'''
		self.check_grid_coordinates(grid_coo)
		return grid_coo[0] * self.n_cols + grid_coo[1]

	def grid_coordinates_to_coordinates(self, grid_coo):
		'''
		Return the physical coordinates of the center of the pixel identified by "grid_coo" (row id, column id)
		'''
		self.check_grid_coordinates(grid_coo)
		return ((self.y_ticks[grid_coo[0]] + self.y_ticks[grid_coo[0] + 1])/2, (self.x_ticks[grid_coo[1]] + self.x_ticks[grid_coo[1] + 1])/2)
	
	def check_grid_coordinates(self, grid_coo):
		'''
		Check if grid coordinates are valid. Raises IndexError if not
		'''
		if grid_coo[0] >= self.n_rows or grid_coo[0] < 0:
			raise IndexError('row id outside the range')
		if grid_coo[1] >= self.n_cols or grid_coo[1] < 0:
			raise IndexError('column id outside the range')

	def check_coordinates(self, coo):
		'''
		Check if physical coordinates are valid. Raises IndexError if not
		'''
		if coo[0] > self.map_size[0] or coo[0] < 0:
			raise IndexError('y coordinate outside the range')
		if coo[1] > self.map_size[1] or coo[1] < 0:
			raise IndexError('x coordinate outside the range')

	def check_node_id(self, node_id):
		'''
		Check if node id is valid. Raises IndexError if not
		'''
		if node_id >= self.N or node_id < 0:
			raise IndexError('Node id outside the range')

	def plot_indegrees(self):
		indegrees = self.A.transpose().dot(np.ones((self.N, )))
		self.plot(indegrees)

	def plot_outdegrees(self):
		outdegrees = self.A.dot(np.ones((self.N, )))
		self.plot(outdegrees)

	def plot(self, values, source=None):
		if len(values) != self.N:
			raise ValueError('Needs a value for each pixel')
		
		z = np.reshape(values, self.shape)
		
		plt.figure()
		plt.imshow(z, cmap='RdYlGn', interpolation='nearest')
		plt.colorbar()
		if source is not None:
			plt.sc
		plt.show()