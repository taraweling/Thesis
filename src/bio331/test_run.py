import run
import utils

toy_nodes = ['v1', 'v2', 'v3', 'v5', 'v6', 'v8', 'v7', 'v4']
toy_edges = [('v1', 'v2'), ('v1', 'v3'), ('v2', 'v3'), ('v3', 'v5'), ('v5', 'v6'), ('v6', 'v8'), ('v8', 'v7'), ('v6', 'v7'), ('v5', 'v7'), ('v5', 'v4')]
toy_adjlist = {'v1': {'v3', 'v2'}, 'v2': {'v1', 'v3'}, 'v3': {'v1', 'v5', 'v2'}, 'v5': {'v7', 'v3', 'v6', 'v4'}, 'v6': {'v8', 'v5', 'v7'}, 'v8': {'v7', 'v6'}, 'v7': {'v8', 'v5', 'v6'}, 'v4': {'v5'}}
toy_paths = run.all_pairs_shortest_paths(toy_adjlist)
toy_B = {('v1', 'v2'): 1.0, ('v1', 'v3'): 6.0, ('v2', 'v3'): 6.0, ('v3', 'v5'): 15.0, ('v5', 'v6'): 7.5, ('v6', 'v8'): 3.5, ('v8', 'v7'): 3.5, ('v6', 'v7'): 1.0, ('v5', 'v7'): 7.5, ('v5', 'v4'): 7.0}
toy_partitions = [[['v1', 'v2', 'v3', 'v5', 'v6', 'v8', 'v7', 'v4']], [['v1', 'v2', 'v3'], ['v4', 'v5', 'v6', 'v7', 'v8']], [['v1', 'v2', 'v3'], ['v5', 'v6', 'v7', 'v8'], ['v4']], [['v1', 'v2', 'v3'], ['v4'], ['v5'], ['v6', 'v7', 'v8']], [['v4'], ['v5'], ['v6', 'v7', 'v8'], ['v1'], ['v2', 'v3']], [['v4'], ['v5'], ['v6', 'v7', 'v8'], ['v1'], ['v2'], ['v3']], [['v4'], ['v5'], ['v1'], ['v2'], ['v3'], ['v8'], ['v6', 'v7']], [['v4'], ['v5'], ['v1'], ['v2'], ['v3'], ['v8'], ['v6'], ['v7']]]
for i in range(len(toy_partitions)):
	for j in range(len(toy_partitions[i])):
		toy_partitions[i][j] = sorted(toy_partitions[i][j])

other_nodes = ['v1', 'v2', 'v3', 'v5', 'v6', 'v8', 'v7', 'v4']
other_edges = [('v1', 'v2'), ('v1', 'v3'), ('v3', 'v5'), ('v5', 'v6'), ('v8', 'v7'), ('v5', 'v4'), ('v2', 'v4'), ('v6', 'v4'), ('v6', 'v7')]
other_adjlist = {'v1': {'v2', 'v3'}, 'v2': {'v1', 'v4'}, 'v3': {'v1', 'v5'}, 'v5': {'v4', 'v3', 'v6'}, 'v6': {'v4', 'v7', 'v5'}, 'v8': {'v7'}, 'v7': {'v8', 'v6'}, 'v4': {'v2', 'v6', 'v5'}}
other_paths = run.all_pairs_shortest_paths(other_adjlist)
other_B = {('v1', 'v2'): 4.5, ('v1', 'v3'): 4.5, ('v3', 'v5'): 7.5, ('v5', 'v6'): 7.5, ('v8', 'v7'): 7.0, ('v5', 'v4'): 3.0, ('v2', 'v4'): 7.5, ('v6', 'v4'): 7.5, ('v6', 'v7'): 12.0}
other_partitions = [[['v1', 'v2', 'v3', 'v5', 'v6', 'v8', 'v7', 'v4']], [['v1', 'v2', 'v3', 'v4', 'v5', 'v6'], ['v7', 'v8']], [['v7', 'v8'], ['v1', 'v2', 'v3'], ['v4', 'v5', 'v6']], [['v7', 'v8'], ['v4', 'v5', 'v6'], ['v1', 'v3'], ['v2']], [['v7', 'v8'], ['v4', 'v5', 'v6'], ['v2'], ['v1'], ['v3']], [['v7', 'v8'], ['v2'], ['v1'], ['v3'], ['v5'], ['v4', 'v6']], [['v2'], ['v1'], ['v3'], ['v5'], ['v4', 'v6'], ['v8'], ['v7']], [['v2'], ['v1'], ['v3'], ['v5'], ['v8'], ['v7'], ['v6'], ['v4']]]
for i in range(len(other_partitions)):
	for j in range(len(other_partitions[i])):
		other_partitions[i][j] = sorted(other_partitions[i][j])

def test_toy_v3_v5():
	yours = run.single_edge_betweenness(toy_nodes,toy_paths,'v3','v5')
	assert yours == 15

def test_toy_v5_v7():
	yours = run.single_edge_betweenness(toy_nodes,toy_paths,'v5','v7')
	assert yours == 7.5

def test_toy_v8_v7():
	yours = run.single_edge_betweenness(toy_nodes,toy_paths,'v8','v7')
	assert yours == 3.5

def test_toy_v7_v8():
	yours = run.single_edge_betweenness(toy_nodes,toy_paths,'v7','v8')
	assert yours == 3.5

def test_other_v3_v5():
	yours = run.single_edge_betweenness(other_nodes,other_paths,'v3','v5')
	assert yours == 7.5

def test_other_v1_v2():
	yours = run.single_edge_betweenness(other_nodes,other_paths,'v1','v2')
	assert yours == 4.5

def test_other_v4_v5():
	yours = run.single_edge_betweenness(other_nodes,other_paths,'v4','v5')
	assert yours == 3

def test_other_v7_v8():
	yours = run.single_edge_betweenness(other_nodes,other_paths,'v7','v8')
	assert yours == 7

def test_toy_B():
	yours = run.edge_betweenness(toy_nodes,toy_edges,toy_paths)
	# the dictionary can have single edges (a,b) or both (a,b) and (b,a).
	assert len(yours) == len(toy_B) or len(yours) == 2*len(toy_B)
	for key in yours.keys():
		e = (key[0],key[1])
		erev = (key[1],key[0])

		assert e in toy_B or erev in toy_B # the edge must be in the dictionary
		if e in toy_B:
			assert toy_B[e] == yours[key]
		else:
			assert toy_B[erev] == yours[key]
	return

def test_other_B():
	yours = run.edge_betweenness(other_nodes,other_edges,other_paths)
	# the dictionary can have single edges (a,b) or both (a,b) and (b,a).
	assert len(yours) == len(other_B) or len(yours) == 2*len(other_B)
	for key in yours.keys():
		e = (key[0],key[1])
		erev = (key[1],key[0])

		assert e in other_B or erev in other_B # the edge must be in the dictionary
		if e in other_B:
			assert other_B[e] == yours[key]
		else:
			assert other_B[erev] == yours[key]
	return

def test_toy_partition_num():
	toy_nodes = ['v1', 'v2', 'v3', 'v5', 'v6', 'v8', 'v7', 'v4']
	toy_edges = [('v1', 'v2'), ('v1', 'v3'), ('v2', 'v3'), ('v3', 'v5'), ('v5', 'v6'), ('v6', 'v8'), ('v8', 'v7'), ('v6', 'v7'), ('v5', 'v7'), ('v5', 'v4')]
	toy_adjlist = {'v1': {'v3', 'v2'}, 'v2': {'v1', 'v3'}, 'v3': {'v1', 'v5', 'v2'}, 'v5': {'v7', 'v3', 'v6', 'v4'}, 'v6': {'v8', 'v5', 'v7'}, 'v8': {'v7', 'v6'}, 'v7': {'v8', 'v5', 'v6'}, 'v4': {'v5'}}
	yours = run.GN(toy_nodes.copy(),toy_edges.copy(), toy_adjlist.copy())

	assert len(yours) == len(toy_partitions)

def test_toy_partition_1():
	assert check_toy_partition(1)
	
def test_toy_partition_2():
	assert check_toy_partition(2)

def check_toy_partition(i):
	toy_nodes = ['v1', 'v2', 'v3', 'v5', 'v6', 'v8', 'v7', 'v4']
	toy_edges = [('v1', 'v2'), ('v1', 'v3'), ('v2', 'v3'), ('v3', 'v5'), ('v5', 'v6'), ('v6', 'v8'), ('v8', 'v7'), ('v6', 'v7'), ('v5', 'v7'), ('v5', 'v4')]
	toy_adjlist = {'v1': {'v3', 'v2'}, 'v2': {'v1', 'v3'}, 'v3': {'v1', 'v5', 'v2'}, 'v5': {'v7', 'v3', 'v6', 'v4'}, 'v6': {'v8', 'v5', 'v7'}, 'v8': {'v7', 'v6'}, 'v7': {'v8', 'v5', 'v6'}, 'v4': {'v5'}}
	yours = run.GN(toy_nodes.copy(),toy_edges.copy(), toy_adjlist.copy())
	
	assert len(yours) > i
	your_p = yours[i]
	sol_p = toy_partitions[i]
	for p in your_p:
		p = sorted(list(p))
		assert p in sol_p
	return True

def test_other_partition_num():
	other_nodes = ['v1', 'v2', 'v3', 'v5', 'v6', 'v8', 'v7', 'v4']
	other_edges = [('v1', 'v2'), ('v1', 'v3'), ('v3', 'v5'), ('v5', 'v6'), ('v8', 'v7'), ('v5', 'v4'), ('v2', 'v4'), ('v6', 'v4'), ('v6', 'v7')]
	other_adjlist = {'v1': {'v2', 'v3'}, 'v2': {'v1', 'v4'}, 'v3': {'v1', 'v5'}, 'v5': {'v4', 'v3', 'v6'}, 'v6': {'v4', 'v7', 'v5'}, 'v8': {'v7'}, 'v7': {'v8', 'v6'}, 'v4': {'v2', 'v6', 'v5'}}
	yours = run.GN(other_nodes.copy(),other_edges.copy(),other_adjlist.copy())
	assert len(yours) == len(other_partitions)

def test_other_partition_1():
	assert check_other_partition(1)

def check_other_partition(i):
	other_nodes = ['v1', 'v2', 'v3', 'v5', 'v6', 'v8', 'v7', 'v4']
	other_edges = [('v1', 'v2'), ('v1', 'v3'), ('v3', 'v5'), ('v5', 'v6'), ('v8', 'v7'), ('v5', 'v4'), ('v2', 'v4'), ('v6', 'v4'), ('v6', 'v7')]
	other_adjlist = {'v1': {'v2', 'v3'}, 'v2': {'v1', 'v4'}, 'v3': {'v1', 'v5'}, 'v5': {'v4', 'v3', 'v6'}, 'v6': {'v4', 'v7', 'v5'}, 'v8': {'v7'}, 'v7': {'v8', 'v6'}, 'v4': {'v2', 'v6', 'v5'}}

	yours = run.GN(other_nodes.copy(),other_edges.copy(), other_adjlist.copy())
	assert len(yours) > i
	your_p = yours[i]
	sol_p = other_partitions[i]
	for p in your_p:
		p = sorted(list(p))
		assert p in sol_p
	return True

	
