def get_center(p1,p2):
    """
    Calculate the center point between two points p1 and p2.
    
    Args:
        p1 (tuple): A tuple representing the first point (x1, y1).
        p2 (tuple): A tuple representing the second point (x2, y2).
    """
    return ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)

def get_distance(c1,c2):
    """
    Calculate the Euclidean distance between two points c1 and c2.
    
    Args:
        c1 (tuple): A tuple representing the first point (x1, y1).
        c2 (tuple): A tuple representing the second point (x2, y2).
    """
    return ((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2) ** 0.5

def save_adj_list_as_txt(adj_list, filename):
	"""
	Save the adjacency list to a text file.
	
	Args:
		adj_list (dict): The adjacency list to save.
		filename (str): The name of the file to save the adjacency list to.
	"""
	with open(filename, 'w') as f:
		for node, edges in adj_list.items():
			for edge in edges:
				edges_str = str(node[0])+ "\t" + str(edge[0]) + "\t" + str(edge[1]) + "\n"
				print(edges_str)
				f.write(edges_str)