from __future__ import division
import numpy as np
import random
import matplotlib.pyplot as plt
from math import sqrt
from scipy.sparse.csgraph import dijkstra
from grid_manager import Grid

class Landmarks:
	
	def __init__(self, grid, n_landmarks=0.01, verbose=True, criterion='uniform'):
		'''
		 - grid is a Grid object representing the landscape

		 - n_landmarks is the number of landmark. If it is bigger than 1, it is rounded to the nearest integer, 
		   if it is in (0, 1), it is interpreted as the fraction of nodes to use as landmarks
		'''
		self.verbose = verbose
		self.G = grid
		self.set_n_landmarks(n_landmarks)
		if criterion=='qualities':
			self.choose_landmarks_from_qualities()
		elif criterion=='uniform':
			self.choose_landmarks()
		else:
			print ("Invalid criterion, choosing landmarks uniformly.")
			self.choose_landmarks()


	def set_n_landmarks(self, n_landmarks):
		'''
		Set self.n_landmarks, the number of landmarks to use, based on the parameter n_landmarks

			If n_landmark it is bigger than 1, it is rounded to the nearest integer, 
			if it is in (0, 1), it is interpreted as the fraction of nodes to use as landmarks
		'''
		N = self.G.N
		if (n_landmarks <= 0 or n_landmarks > N):
			raise ValueError('n_landmarks must be in the interval (0, number of nodes]')
		elif n_landmarks >= 1:
			self.n_landmarks = int(round(n_landmarks))
		else:
			self.n_landmarks = int(round(N * n_landmarks))

	def choose_landmarks(self):
		'''
		Select landmarks uniformly distributed on the grid
		'''

		n_rows, n_cols = self.G.shape
		N = self.G.N

		l_cols = int(max(1, min(self.n_landmarks, round(sqrt(n_cols*self.n_landmarks/n_rows)))))
		
		l_full_rows = int(self.n_landmarks//l_cols)
		n_extra = int(self.n_landmarks - l_cols*l_full_rows)
		
		l_rows = l_full_rows
		if n_extra > 0:
			l_rows += 1
			print('Warning: Landmarks are not uniformly spread. Consider using {} or {} landmarks instead of {}'.format(self.n_landmarks - n_extra, self.n_landmarks - n_extra + l_cols, self.n_landmarks))
		
		rows_indices = [int((2*i+1) * n_rows / l_rows / 2) for i in range(l_rows)]
		col_indices = [int((2*i+1) * n_cols / l_cols / 2) for i in range(l_cols)]
		extra_indices = [int((2*i+1) * n_cols / n_extra / 2) for i in range(n_extra)]
		
		self.landmarks = np.zeros((self.n_landmarks,), dtype=np.int32)
		i = 0
		for r in rows_indices[:l_full_rows]:
			for c in col_indices:
				self.landmarks[i] = self.G.grid_coordinates_to_node_id((r,c))
				i += 1
		for c in extra_indices:
			self.landmarks[i] = self.G.grid_coordinates_to_node_id((rows_indices[-1],c))
			i += 1
	
	def choose_landmarks_from_qualities(self):
		 highest_quality_nodes = np.argsort(self.G.qualities)[-self.n_landmarks*100:-1][::-1]
		 self.landmarks = random.sample(highest_quality_nodes, self.n_landmarks)

	def plot_landmarks(self):

		values = np.zeros((self.G.N,))
		values[self.landmarks] = 1
		self.G.plot(values)

	def plot_landmark(self, landmark_id):

		values = np.zeros((self.G.N,))
		print(landmark_id)
		print(self.landmarks[landmark_id])
		values[self.landmarks[landmark_id]] = 1
		self.G.plot(values)

	def plot_landmark_similarities(self, landmark_id):
		
		affinities = self.similarities_L2all[landmark_id, :]
		self.G.plot(affinities)

	def similarities_to_landmarks(self, min_affinity=0, affinity_to_cost=None, distance_transformation=None, distance_to_similarity=None):
		'''
		Compute the affinities from landmarks to all nodes, and from all nodes to landmarks
		Store it in self.L2all (horizontal: n_landmarks x n_nodes) and self.all2L (vertical: n_nodes x n_landmarks)
		
		N.B. : To avoid confusion, cost and affinities are properties of edges, while distances and similarities are properties of pair of nodes that are not necessarily directly connected.
		
		Parameters
		----------
			min_affinity: float (default = 1e-3)
			Affinities below this value will not be considered as possible steps.
			They will be associated with an infinite cost, will have an infinite distance to all landmark, 
			and a 0 similarity to all node.

			affinity_to_cost: function handle (default = None)
			Function that defines the transformation of affinities into costs. 
			It must accept a numpy vector of floats (the affinities) and return a vector of the same size (the costs)
			if None, cost are defined as 1 / affinities.

			distance_transformation: function handle (default = None)
			Function that would be applied to the results the dijkstra distances computation.
			self.distances_L2all = distance_transformation(dijkstra_results).
			If distance_transformation is None, self.distances_L2all = dijkstra_results.
			The same goes for self.distances_all2L.
			distance_transformation=np.log can be a good choice

			distance_to_similarity: function handle (default = None)
			Function that defines the transformation of distances into similarities.
			if None, similarities are defined as (max_distance - distance) / max_distance.

		'''

		def finite_max(a):
			''' Returns the maximum of an array, excluding infinte values
			'''
			return np.amax(a[np.isfinite(a)])
		
		
		A = self.G.A.copy() # make a copy of the affinity matrix

		# Transformation of affinities to cost
		impossible_steps = A.data < min_affinity
		if affinity_to_cost == None: 	# default affinity to cost transformation
			A.data[:] = 1. / A.data
		else: 							# user specified affinity to cost transformation
			A.data[:] = affinity_to_cost(A.data)
		A.data[impossible_steps] = np.inf

		
		# From landmarks to all points
		##############################
		if self.verbose:
			print("Compute similarities from landmarks to all nodes...")

		# Shortest-path computation with Dijkstra
		self.distances_L2all = dijkstra(A, indices=self.landmarks)

		# Optional distance transformation
		if distance_transformation != None: 
			self.distances_L2all = distance_transformation(self.distances_L2all)

		# Distance to similarity transformation
		if distance_to_similarity == None: 	# default distance to similarity transformation
			self.similarities_L2all = (finite_max(self.distances_L2all) - self.distances_L2all)/finite_max(self.distances_L2all)
			self.similarities_L2all[self.distances_L2all == np.inf] = 0
		else:								# user specified distance to similarity transformation
			self.similarities_L2all = distance_to_similarity(self.distances_L2all)

		# From all points to landmarks
		##############################
		if self.verbose:
			print("Compute similarities from all nodes to landmarks...")

		# Shortest-path computation with Dijkstra
		self.distances_all2L = np.transpose(dijkstra(A.transpose(), indices=self.landmarks))

		# Optional distance transformation
		if distance_transformation != None: 
			self.distances_all2L = distance_transformation(self.distances_all2L)

		# Distance to similarity transformation
		if distance_to_similarity == None: 	# default distance to similarity transformation
			self.similarities_all2L = (finite_max(self.distances_all2L) - self.distances_all2L)/finite_max(self.distances_all2L)
			self.similarities_all2L[self.distances_all2L == np.inf] = 0
		else:								# user specified distance to similarity transformation
			self.similarities_all2L = distance_to_similarity(self.distances_all2L)

	def similarities_from(self, sources):
		'''
		Compute the similarities based on landmarks, from the sources to every nodes in the grid.
		sources is either a single int or a list of int, representing the indices of the sources.
		'''

		return self.similarities_all2L[sources, :].dot(self.similarities_L2all)

	def similarities_to(self, destinations):
		'''
		Compute the similarities based on landmarks, from every nodes to the destinations.
		destinations is either a single int or a list of int, representing the indices of the destinations.
		'''

		return self.similarities_all2L.dot(self.similarities_L2all[:, destinations])

	def habitat_functionalities(self):
		'''
		Compute habitat functionalities of all the pixels as sums of similarities from
		an origin pixel to all destination pixels, weighted by the destinations' quality values.
		'''
		L2all_Q = self.similarities_L2all.dot(self.G.qualities)
		hf = self.similarities_all2L.dot(L2all_Q)
		hf[hf == 0] = np.nan
		return hf

	def harmonic_approximation(self):
		all2L_Q = self.similarities_all2L.dot(self.G.qualities[self.landmarks])
		all2L_Q[all2L_Q==0] = np.nan
		return all2L_Q.tolist()

	#  def closeness_centrality(self, weight=None):
	# 	'''
	# 	Compute the closeness centrality (approximated with landmarks)
	# 	This centrality is based on shortest path distances, and assumes that self.distances_all2L and self.distances_L2all contains the shortest path distances from and to all landmarks.
	# 	The exact 
	# 	'''