[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/J43xVlF8)
# Community Detection on Food Webs

**Purpose:** The goal of this assignment is to implement the Girvan-Newman algorithm for community detection and apply it to the Chesapeake Bay Food Web. This assignment builds directly on the edge betweenness code you wrote in Lab 5. It is the last programming assignment before the research project, so it is less structured than previous assignments and labs to give you practice designing your own algorithms.

**Criteria for Success:** There are six tasks in this assignment (Tasks A-F); Task C is the biggest task, followed by Tasks E & F that ask you to visualize and analyze a food web. Some tasks have tests, and all tasks have points assigned to them (some will be manually graded). **There are 16 total points: 0-8 is Missing, 9-11 is Partial, 12-14 is Complete, and 15-16 is Excellent.** 

| Task | Pytests | Points | Notes | 
| -- | -- | -- | -- | 
| Task A: Read an Edge File | 0 | 1 | Previous code for reading edge info | 
| Task B1: Single Edge Betweenness | 8 | 1 | Lab 5 |
| Task B2: Edge Betweenness | 2 | 2 | One point for each `pytest` test |
| Task C: Girvan-Newman | 5 | 5 | One point for each step |
| Task D: Viz Example Network | 0 | 1 | Function in `utils.py` |
| Task E: Viz Bay Network | 0 | 2 | P2 for reading node info |
| Task F: Compare Clusterings | 0 | 2 | Write a comment in `run.py` |
| Readable Code | 0 | 2 | |
| **Total** | **15** | **16** | 

The "Readable Code" criterion means that your code is organized in a coherent way and you use comments throughout your code to explain any tricky parts. Specifically, I am looking for:
- Consistent conventions and clear variable names (however you choose to comment, be consistent).
- Use of whitespace and linebreaks to separate sections of code within a function.
- In-line comments that explain the rationale of any tricky parts. 
- Removal of comment stubs that no longer make sense after you have completed the assignment.
- [Attribution of any collaborators](https://reed-compbio-classes.github.io/bio331-syllabus/doc/policies/#collaboration-policy), [websites](https://reed-compbio-classes.github.io/bio331-syllabus/doc/policies/#online-resources--generative-ai-policy), or [generative AI](https://reed-compbio-classes.github.io/bio331-syllabus/doc/policies/#online-resources--generative-ai-policy) (include the prompts in this case).

At any point you can type `pytest` in the terminal to see which tests pass. These tests are located in `test_run.py`, so you can see which examples pass and fail and debug your program.

:exclamation: Do not modify `test_run.py`. If you’d like to experiment with your own tests, create a new file beginning with `test_`, which pytest will also detect.

**Resources:** Refer to the [Bio131 Python Crashcourse](https://reed-compbio-classes.github.io/python-crashcourse/) for Python syntax refreshers. Take a look at the Files and the Manipulating Strings sections. [PythonTutor](http://pythontutor.com/) can also help you with interactive debugging. 

**Collaboration is encouraged!** You may also copy code from previous assignments and labs (as long as you cite where the code is from in your comments). 

:alien: I used ChatGPT to identify and clarify ambiguous instructions in this assignment. 

:alien: You can ask generative AI simple Python syntax questions, but do not use it for substantial coding questions (such as writing entire functions). You must include the prompts for any ChatGPT-generated queries as comments in your code. [See the programming policy](https://reed-compbio-classes.github.io/bio331-syllabus/doc/policies/#programming-assignments) on the syllabus. 

## The Goal

Your goal is to implement the Girvan-Newman algorithm and visualize the Chesapeake Bay Food Web with different numbers of communities. The Girvan-Newman algorithm is from the paper [Community Structure in Social and Biological Networks](https://arxiv.org/pdf/cond-mat/0112110v1.pdf) (PNAS 2002).    

**For this assignment, all graphs are undirected, unweighted, and connected.**

This divisive algorithm has the following steps:

1. Assign all nodes to a single cluster.
2. Calculate the _edge betweenness_ of all edges in the network.  
3. Remove the edge with the highest betweenness.
4. If removing the edge `(u,v)` means that `u` and `v` are no longer in the same connected component, replace the single cluster containing `u` and `v` with two clusters, one that contains `u` and one that contains `v`.
5. Repeat from Step 2 until no edges remain to remove.

The figure below shows an example graph with eight nodes and the first two iterations of Girvan-Newman that remove ($v_3$,$v_5$) and ($v_4$,$v_5$), respectively.

![GN](figs/GN.png)

Your method should output a list of partitions, where each new partition increases the number of clusters by one.

:question: Have a question at any point? Ask Anna!  This assignment has many components, and there's a good chance that Anna can provide additional clarifications to the class.

:construction: Functions labeled with :construction: are written for you in `utils.py`. You can use these or you can choose to implement them on your own, if you want a challenge.

## Preliminaries: Run `run.py`

**Activate your virtual environment:** In VSCode, open the Command Palette (`Ctrl/Cmd` + `Shift` + `P`) and select your Python interpreter for your environment (`~/bio331-venv/bin/activate` for venv or the `bio331` conda environment). You must select your environment each time you open a new lab or programming assignment in VSCode.

**Change to the directory that contains `run.py`:** When you open the terminal in VSCode, it might not be in your local repository location. Change directories (with the `cd` command) until you are in your repository (and you can see `run.py` when you type `ls`, which means "list all files & folders").

**Run the provided code:** Run the provided `run.py` file. It should produce no output but also no errors. The `run.py` file uses helper functions from `utils.py`. For this assignment, you will write new code in `run.py` but you will copy code into both `run.py` and `utils.py`.

**One Handy Trick: `sys.exit()`** There are lots of iterations and loops here; sometimes when debugging you want to print something within a loop and then quit right away. You can type `sys.exit()` on any line and Python will immediately quit. For example, the code snippet below will quit when `i > 10`, even if it's in an infinite loop:

```
i = 0
while True: # oops this is an infinite loop!
  i+=1
  if i > 10:
    sys.exit()
```

## :star: **Task A**: Read the Edge File

**For Tasks A-C, uncomment the lines in `main()` that call the relevant functions as you work.**

There are two graphs included in P4. The first graph, `files/example-edges.txt`, contains the edges from the example in the image shown above; start with this one. The second graph is the Chesapeake Food Web, and you will work on that in Task D.

For Task A, write a `read_edges()` function to read `files/example-edges.txt`. The function takes a filename (a string). Return (1) a list of nodes, (2) a list of 2-element lists for edges, and (3) a dictionary-based adjacency list. At this point you have seen code to read a graph many times; be sure to cite where you copied code from. This is a connected graph, so all nodes appear in the edgelist.

## :star: **Task B**: Implement Edge Betweenness 

### B1. Calculate the edge betweenness for _one_ edge.  

You've done this already! Copy your code from Lab 5 (or the solution), remember to cite in the comments where the code is from. Also, you will need to copy `shortest_paths_ties()` and `get_paths()` from Lab 5's `utils.py` file; put them in P4's `utils.py` file so you don't need to update code from Lab 5. 

**There are 8 pytests for the `single_edge_betweenness()` function.** Four of the pytests are from the example graph above and four are from another graph (called "other" in the tests):

| Example Graph (used in instructions) | Test Graph (used for additional pytests) | 
| -- | -- |
| ![example-graph](figs/example-graph.png) | ![test-graph](figs/test-graph.png) |

### B2. Calculate edge betweenness for all edges in the graph.

Now, write an `edge_betweenness()` function that takes a node list/set, an edges list, and the `paths` variable from `utils.all_pairs_shortest_paths(adjlist)`. It returns a dictionary where the keys are edges (represented as _tuples_) and the values are the edge betweenness centrality (floats). Within this function, you will call the `single_edge_betweenness()` function from Task B1.

The final betweenness dictionary will have _(e,B(e))_ key-value pairs. To store edges as _keys_ in a dictionary, use a tuple (e.g., `(u,v)` or `('v1','v2')`) rather than 2-element lists. (`[u,v]`). Note the parentheses instead of lists. Tuples can be indexed into just like lists:

```
edge = ('a','b')
print(edge[0])
print(edge[1])
```

The difference is that you cannot change the contents of a tuple, so they can be hashed as keys. (This is a minor detail, but it's helpful to know about if you're interested in this kind of thing).

**Expected Output:** In the first iteration of the algorithm, the edge betweenness dictionary contains the following values. Note that the edge is represented only once - the pytest tests will check if exactly one of each edge is present (either `('v1', 'v2')` or `('v2', 'v1')` should be in your dictionary, but not both). The order of the keys does not matter because dictionaries are unordered.

```
key        value
('v1', 'v2') 1.0
('v1', 'v3') 6.0
('v2', 'v3') 6.0
('v3', 'v5') 15.0
('v5', 'v6') 7.5
('v5', 'v7') 7.5
('v5', 'v4') 7.0
('v6', 'v8') 3.5
('v6', 'v7') 1.0
('v8', 'v7') 3.5
```

**There are two pytests for the `edge_betweenness()` function.** One pytest checks that your betweenness dictionary contains the values above, and the other pytest checks that your dictionary calculates edge betweenness for the dark blue "other" graph.

## :star: **Task C**: Implement the Girvan-Newman Algorithm

Now, write a `GN()` function that takes a node list/set, an edge list, and an adjacency list and returns a list of lists that describes all of the clusters at each step of the algorithm. Let's revisit the steps for the Girvan-Newman Algorithm:

### C1. Assign all nodes to a single cluster.

Create a `partitions` variable, which will be a list of lists.  Each element in `partitions` will be _one_ grouping.  This starts with all nodes in a single group. Initialize the `partitions` variable with the `nodes` list.

```
partitions = [[nodes]]
```

When you print `partitions` you will see a nested list like this:

```
[[['v1', 'v2', 'v3', 'v5', 'v6', 'v8', 'v7', 'v4']]]
```

This is expected! It means that you start with exactly one clustering, `[['v1', 'v2', 'v3', 'v5', 'v6', 'v8', 'v7', 'v4']]`, and that contains exactly one cluster, `['v1', 'v2', 'v3', 'v5', 'v6', 'v8', 'v7', 'v4']`. Your `GN()` function will **return** the `partitions` variable after you've appended other clusterings to it. An example of all partitions appear later in this section.

### C2. Calculate the _edge betweenness_ of all edges in the network.  

You've done it! Call the code you wrote for Task B.

### C3. Remove the edge with the highest betweenness

You should get `('v3','v5')` in the first iteration; you need to remove this edge from the graph. Note that the edge is a _tuple_ in the edge betweenness dictionary (from B2), but you need to convert it back to a list (e.g., `['v3','v5']`) when removing from the edge list and adjlist.

:construction: There are two functions in `utils.py` that will come in handy.  The first removes an edge from an edgelist, the second removes an edge from an adjacency list. Both of these functions change the data structures _in-place_, meaning that you don't need to assign them to new variables. Try running these examples (either in the `main()` function or in the `GN()` function) to see what happens. Comment them out when you understand how they work.

```
edgelist = [['a','b'],['b','c'],['d','a']]
utils.remove_from_edgelist(['a','b'],edgelist)
print(edgelist)
utils.remove_from_edgelist(['a','d'],edgelist)
print(edgelist)


adjlist = {'a':['b','c'],'b':['a','c','d'],'c':['a','b','c'],'d':['b']}
utils.remove_from_adjlist(['a','b'],adjlist)
print(adjlist)
utils.remove_from_adjlist(['d','b'],adjlist)
print(adjlist)
```

:bulb: You need to remove the edge from _both_ the edge list and the adjacency list - keep these two data structures consistent.

### C4. If removing an edge divided a group, make a new partition.

If removing the edge `(u,v)` means that `u` and `v` are no longer in the same connected component, replace the single cluster containing `u` and `v` with two clusters, one that contains `u` and one that contains `v`. There is some code in `utils.py` to help you here - to find the connected component containing a node and the split a partition. Your task is to check to see if `u` and `v` are in the same component or not, and if they are not then split the partition. 

:construction: The `utils.conncomp()` function that takes an adjacency list and a single node `u` and returns a sorted list of nodes that are in the same connected component as `u`. (It is modified from L4's `check_cycle()` function, in fact). 

:construction: The `utils.split_partition()` function takes the _previous partition_ and the two splits (two lists or sets of nodes). The function returns the updated partition. The previous partition can be accessed using `partitions[len(partitions)-1]`, or even cooler, `partitions[-1]`.  Here's some example code to run (in `main()` or `GN()`; comment it out when done):

```
prev_partition = [['a','b','c','d'],['e','f']] 
split1 = ['a','c']
split2 = ['b','d']
next_partition = utils.split_partition(prev_partition,split1,split2)
print(next_partition)
```

`next_partition` will then be:

```
[['e', 'f'], ['a', 'c'], ['b', 'd']]
```

### C5. Repeat from Step 2 until no edges remain to remove.

This will require a `while` loop, since you don't know when removing an edge results in a split. Your `while` loop will contain steps C2-C4.

:question: What is the terminating condition? You will need to figure this out.

### Return all partitions.

Return a list of partitions, where each partition is itself a list of node groups (lists). Note that `partitions` is a _single_ list, where each element of the list is a partition (or grouping) of _all_ the nodes.  One final partition might be the following (line breaks help visually see the partitions, your parititon might be different):

```
[
[['v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v8']],
[['v1', 'v2', 'v3'], ['v4', 'v5', 'v6', 'v7', 'v8']],
[['v1', 'v2', 'v3'], ['v4'], ['v5', 'v6', 'v7', 'v8']],
[['v1', 'v2', 'v3'], ['v4'], ['v5'], ['v6', 'v7', 'v8']], # group at index 0 is separated in any order
[['v1', 'v2'], ['v3'], ['v4'], ['v5'], ['v6', 'v7', 'v8']],
[['v1'], ['v2'], ['v3'], ['v4'], ['v5'], ['v6', 'v7', 'v8']], # group at index 5 is separated in any order
[['v1'], ['v2'], ['v3'], ['v4'], ['v5'], ['v6'], ['v7', 'v8']],
[['v1'], ['v2'], ['v3'], ['v4'], ['v5'], ['v6'], ['v7'], ['v8']]
]
```

**There are 5 pytests for the `GN()` function.** Two pytests check that the length of `partitions` is correct. There should be one partition for each number of nodes in the graph (since each iteration removes one edge and eventually isolates every node). Two pytests check the first two partitions for the example above and a last pytest checks the first partition in the "other" example. These are unique partitions (as in, there are no ties), so your code should return the same partitions.

## :star: **Task D**: Post the Example Graph Colored by Communities

**For Tasks D & E, you must add lines to `main()` that _call_ the visualization functions.**

All right! You can now make a graph that colors the nodes by partition. Again, there is a `utils.py` function to help you. This function is just like previous functions, but it also takes a partition of nodes. If your `partitions` variable looks like the one above, then you can get a single partition by taking an element from `partitions` (e.g., `partitions[3]` will return the partition of 4 groups). If you have the `partitions` variable and the `edgelist` variable, you can call:

```
utils.viz_example(nodes, edgelist, partitions[3], 'example-4groups.html')
```

![example-communities](figs/example-communities.png)

:exclamation: Ack! You might not have any edges appear!  This happens because the edge list and adjacency list were modified during Girvan–Newman. Since the `edgelist` and `adjlist` variables were modified in-place, they will be the leftovers after Girvan-Newman is done. Before calling `utils.viz_example()`, simply re-read the edges file by calling the `read_edges()` function again to get new versions of the edge list and adjacency list. (There are other ways to handle this too).

:bulb: The colors are randomly generated each time you call this function. Re-run if the color palette is hard to decipher. You can also go into the `viz_example()` function in `utils.py` and hard-code at least 5 colors.

## :star: **Task E**: Apply and Visualize Girvan-Newman on the Chesapeake Bay Food Web

In the last two tasks, you'll add a bit more code to run `GN()` on the Chesapeake Bay food web. You will visualize the communities, _and_ some information about each node. There are two files related to the Chesapeake Bay foodweb:

- `files/foodweb-edges.txt`: list of edges, where taxa have IDs.
- `files/foodweb-nodes.txt`: a three-column file of (a) taxon ID, (b) taxon name, and (c) taxon classification (Benthic, Pelagic, or Unknown).

The node file will be helpful for annotating your graph with node names and node shapes according to their classification.  

You are now ready to run the code on the Chesapeake Bay Food Web.
1. Read the foodweb edges and the foodweb nodes. :construction: P2 may be helpful to read node information.
2. Run Girvan-Newman to generate all partitions.
3. Make a new `viz_foodweb()` function by copying `viz_example()` from the `utils.py` file and modifying it. The function can be in `utils.py` or `run.py`. In `viz_foodweb()`, **annotate the node label** according to the taxon name and the **annotated the node shape** according to its classification: 'Benthic','Pelagic', or 'Unknown'. This information comes from `files/foodweb-nodes.txt`.

:question: What node shapes are available? The [pyvis tutorial](https://pyvis.readthedocs.io/en/latest/tutorial.html) points to a list of shapes from [visjs documentation](https://visjs.github.io/vis-network/docs/network/nodes.html). The shapes with the label inside of it are: ellipse, circle, database, box, text. The ones with the label outside of it are: image, circularImage, diamond, dot, star, triangle, triangleDown, hexagon, square and icon. You specify node shape or node label when you call the `add_node()` function (e.g., shape='diamond').

:question: Note that this function calls other functions in the `utils.py` file. In those cases, you need to change `rgb_to_hex()` to `utils.rgb_to_hex()` so Python knows where to look for the function. 

4. Generate two HTML files for visualization: one showing two clusters and another showing four clusters. Give them different filenames that correspond to the number of clusters (e.g., `foodweb-2clusters.html`).

## :star: **Task F**: Compare the Unweighted Clustering to the Paper Result

Your graph was generated using an unweighted graph that captures the highest-weight edges of the original Chesapeake Bay Food Web. Running the weighted version produced this result from the [original paper by Girvan & Newman](https://www.pnas.org/doi/10.1073/pnas.122653799). Here's an excerpt from their paper:

_We have also applied our algorithm to a food web of marine organisms living in the Chesapeake Bay, a large estuary on the east coast of the United States. This network was originally compiled by Baird and Ulanowicz ([27](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.2307/1943071)) and contains 33 vertices representing the ecosystem's most prominent taxa. Most taxa are represented at the species or genus level, although some vertices represent larger groups of related species. Edges between taxa indicate trophic relationships—one taxon feeding on another. Although relationships of this kind are inherently directed, we here ignore direction and consider the network to be undirected._

_Applying our algorithm to this network, we find two well defined communities of roughly equal size, plus a small number of vertices that belong to neither community (see Fig. 7). As Fig. 7 shows, the split between the two large communities corresponds quite closely with the division between pelagic organisms (those that dwell principally near the surface or in the middle depths of the bay) and benthic organisms (those that dwell near the bottom). Interestingly, the algorithm includes within each group organisms from a variety of different trophic levels. This finding contrasts with other techniques that have been used to analyze food webs ([28](https://www.jstor.org/stable/3546748)), which tend to cluster taxa according to trophic level rather than habitat. Our results seem to imply that pelagic and benthic organisms in the Chesapeake Bay can be separated into reasonably self-contained ecological subsystems. The separation is not perfect: a small number of benthic organisms find their way into the pelagic community, presumably indicating that these species play a substantial role in the food chains of their surface-dwelling colleagues. This finding suggests that the simple traditional division of taxa into pelagic or benthic may not be an ideal classification in this case._

![tree](figs/foodweb-tree.png)

In a comment near the end of `run.py`, write a short paragraph about the similarities and differences you observe in your network compared to the published result.

## Optional Challenge: Use the weighted food web!

There is one more file for you to work with - `files/foodweb-edges-weighted.txt` contains three columns: node1, node2, and weight. The `files/foodweb-nodes-weighted.txt` includes a few more nodes that had low-weight edges that were removed in the unweighted version.

In order to run your code on a weighted network, you need to do the following two things:

1. Convert the edge weights, where higher is better, to edge _costs_. Typically, this is done by taking the negative log of the weight. You can add `import math` to the top of your file and then calculate the log of a number `x` by calling `math.log(x)`.

![neglog](figs/neglog.png)

2. Implement Dijkstra's algorithm for calculating shortest paths on weighted networks. See the slides from Week 2 to review how to implement Dijkstra's algorithm. :bulb: It is quite similar to the shortest_paths code in the utils file.

3. Add new `weighted_GN()` and `weighted_all_pairs_shortest_paths()` functions to call Dijkstra's instead of the unweighted shortest paths code. 

## Submitting

:star2: **You're Done with Tasks A-F!** Commit and push your changes, then confirm that your updated files appear in your repository in the [Reed Compbio Classes GitHub organization](https://github.com/Reed-Compbio-Classes/).