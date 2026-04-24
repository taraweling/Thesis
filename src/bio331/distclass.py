# code to produce adjacency list of distances from an image
## by miles
from useful import get_distance
import csv

def process_data(filename:str,threshold:float): 
	"""
	Process the data from the given filename.
	
	Args:
		filename (str): The path to the data file.
	"""
	processed_data = [] # List to hold processed data
	adj = {}
	with open(filename, 'r') as file:
		reader = csv.reader(file)
		for node,row in enumerate(reader):
			if ("X" == row[1]): # if the first column is X, skip the row
				continue
			x = row[1]
			y = row[2] ##this is gross and wrong
			processed_data.append([(float(x),float(y)),node]) ## this is just getting the center of each pair of points, 
			# I am not sure that this is the way to do it given the data, will need to be adjusted
		
		for u in processed_data: #for node pair in processed data
			for v in processed_data:
				if u != v: ## not trying to connect a node to itself
					n1,n2 = u[1],v[1] ## nodes get the node numbers
					dist = get_distance(u[0],v[0]) ## get the distance between the two nodes
					if dist <= threshold: #threshold is the distance threshold
						adj.setdefault(n1, []) # if n1 is not in adj, add it
						adj.setdefault(n2, []) # if n2 is not in adj, add it
						edge_uv = [n2,dist] # edge from u to v with distance
						edge_vu = [n1,dist] # edge from v to u with distance
						if edge_vu not in adj[n2]: ## n2 is v
							adj[n2].append(edge_vu) # add edge from v to u
						if edge_uv not in adj[n1]: ## n1 is u
							adj[n1].append(edge_uv)  # add edge from u to v
		return adj #adjacency list of the graph

