"""
Python module for Numerical Integration
29th Summer Institute in the Natural Science and Mathematics
27 April 2023

Gilbert Peralta
Department of Mathematics and Computer Science
College of Science
University of the Philippines Baguio
Email: grperalta@up.edu.ph
Website: https://dmcsweb.upb.edu.ph/~grperalta/
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.patches import Rectangle

#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

UPMaroon = [123/255, 17/255, 19/255]
UPGreen = [1/255, 68/255, 33/255]

class RiemannSum:
	"""
	Class for Riemman sums
	
	Attributes
	----------
		function : callable
			The function to be integrated. Default is f(x) = x.
		left_endpoint : float
			The lower limit of integration. Default value is 0.
		right_endpoint : float
			The upper limit of integration. Default value is 1.
		num_points : int
			The number of points in the partition. Defaut value is 2.
		sample_point : str
			The type of sample points, either "mid", "left", "right", or "random".
			Default is "mid".
		partition : array
			The collection of points subdividing the domain including the endpoints.
			Default is the parition consisting of only the endpoints. To get a uniform
			partition, call self.set_uniform_partition(). To get a random partition, call
			self.set_random_partition().
	"""
	def __init__(self):
		"""
		Class initialization.
		"""
		self.function = lambda x : x
		self.left_endpoint = 0.
		self.right_endpoint = 1.
		self.num_points = 2
		self.sample_point = "mid"
		self.set_uniform_partition()
		
	def set_uniform_partition(self):
		"""
		Returns an array of points corresponding to a uniform partition of the domain of 
		integration based on the number of points (self.num_points).
		"""
		self.partition \
			= np.linspace(self.left_endpoint, self.right_endpoint, self.num_points)	
		self.partition_name = "Uniform"
			
	def set_random_partition(self):
		"""
		Returns an array of points corresponding to a random partition of the domain of 
		integration based on the number of points (self.num_points).
		"""
		self.partition = np.zeros(self.num_points)
		self.partition = np.sort(np.random.rand(self.num_points)*(
			self.right_endpoint - self.left_endpoint) + self.left_endpoint)
		self.partition[0] = self.left_endpoint
		self.partition[-1] = self.right_endpoint 
		self.partition_name = "Random"
		
	def set_sample_points(self):
		"""
		Returns an array of sample points based on the partition (self.partition) and
		the type of sample point (self.sample_point).
		"""
		if self.sample_point == "mid":
			self.sample_points \
				= np.array([0.5*(self.partition[k] + self.partition[k+1]) 
				for k in range(self.num_points-1)])
		if self.sample_point == "left":
			self.sample_points \
				= np.array([self.partition[k] for k in range(self.num_points-1)])
		if self.sample_point == "right":
			self.sample_points \
				= np.array([self.partition[k+1] for k in range(self.num_points-1)])
		if self.sample_point == "random":
			self.sample_points \
				= np.array([np.random.rand()*(self.partition[k+1] - self.partition[k]) 
				+ self.partition[k] for k in range(self.num_points-1)])
		if self.sample_point not in ["mid", "left", "right", "random"]:
			print("UserWarning: Incorrect keyword. Must be either of the following: " 
			+ "'mid', 'left', 'right', 'random'. Performing 'mid' keyword instead.")
			self.sample_point = "mid"
			self.sample_points \
				= np.array([0.5*(self.partition[k] + self.partition[k+1]) 
				for k in range(self.num_points-1)])
		
	def riemann_sum(self):
		"""
		Returns the Riemann sum.
		"""
		self.set_sample_points()
		self.heights = self.function(self.sample_points)
		self.widths = np.array([self.partition[k+1] - self.partition[k] 
			for k in range(self.num_points-1)])
		self.sum = np.dot(self.widths, self.heights)
		return self.sum
		
	def trap_sum(self):
		"""
		Returns the trapezoidal sum.
		"""
		widths = np.array([self.partition[k+1] - self.partition[k] 
			for k in range(self.num_points-1)])
		heights = np.array([0.5*(self.function(self.partition[k+1]) + 
			self.function(self.partition[k])) for k in range(self.num_points-1)])
		_trap_sum = np.dot(widths, heights)
		return _trap_sum
		
	def simp_sum(self):
		"""
		Returns the Simpson sum.
		"""
		widths = np.array([self.partition[k+1] - self.partition[k] 
			for k in range(self.num_points-1)]) / 2
		midpts = np.array([(self.partition[k+1] + self.partition[k]) / 2
			for k in range(self.num_points-1)])
		heights = np.array([(self.function(self.partition[k+1]) + 
			4*self.function(midpts[k]) + self.function(self.partition[k])) / 3 for k in range(self.num_points-1)])
		_simp_sum = np.dot(widths, heights)
		return _simp_sum
		
	def plot(self, gridsize=100, graphcolor=UPMaroon, edgecolor="black", 
		facecolor=UPGreen, alpha=0.5, showheight=True, showgrid=True,
		equalaxis=False):
		"""
		Plots the Riemman sum.
		"""
		if hasattr(self, "sum"):
			pass
		else:
			self.riemann_sum()
		grid = np.linspace(self.left_endpoint, self.right_endpoint, gridsize)
		fig, ax = plt.subplots(figsize=(10, 8))
		for k in range(len(self.partition)-1):
			if showheight:
				ax.plot([self.sample_points[k], self.sample_points[k]], 
					[0, self.function(self.sample_points[k])],
					color="black", lw=0.5, ls="-.")
			ax.add_patch(Rectangle((self.partition[k], 0), 
				self.widths[k], self.heights[k], 
				edgecolor=edgecolor, facecolor=facecolor, alpha=alpha))	
		ax.plot(grid, self.function(grid), color=graphcolor)
		if showgrid:
			ax.grid(which="both", color="lightgray", ls="--")
		ax.axhline(y=0, color="black", lw=0.5)
		ax.axvline(x=0, color="black", lw=0.5)
		ax_title = "Riemann Sum Using {} Partition with {} as Samples: \n {}".format(
			self.partition_name, self._sample_point_name(), self.sum)
		ax.set_title(ax_title)
		if equalaxis:
			ax.axis("equal")
		plt.show()
		
	def norm_partition(self):
		"""
		Returns the norm of the partition (self.partition), that is, the largest
		length of the subinterval induced by the partition.
		"""
		return max([self.partition[k+1] - self.partition[k] 
			for k in range(self.num_points-1)])
		
	def _sample_point_name(self):
		"""
		Returns the following "Left Endpoints", "Right Endpoints", "Midpoints",
		"Random Points" if value of self.sample_point is either "left", "right",
		"mid", "random", respectively.
		"""
		if self.sample_point == "left":
			return "Left Endpoints"
		if self.sample_point == "right":
			return "Right Endpoints"
		if self.sample_point == "mid":
			return "Midpoints"
		if self.sample_point == "random":
			return "Random Points"
			
def logplot(xdata, ydata):
	"""
	Plots the xdata versus ydata in logratihmic scale.
	"""
	fig, ax = plt.subplots(figsize=(10, 8))
	ax.loglog(xdata, ydata, color=UPMaroon, ls="--", marker=".", lw=1.5, ms=15)
	ax.grid(which="major", color="lightgray", ls="--")
	plt.show()