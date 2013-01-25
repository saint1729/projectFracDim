import networkx as nx
import math
G=nx.Graph()
G_=nx.Graph()
colors = [1, 2, 3, 4,  5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];

filepath = raw_input("Enter Filepath: ")
f = open(filepath)
a = [int(n) for n in f.read().split()]
G.add_nodes_from(a)
G_.add_nodes_from(a)

i = 0
b = []
while i<len(a):
    b.append((a[i],a[i+1]))
    b.append((a[i+1],a[i]))
    i = i+2

G.add_edges_from(b)
nodes_colors={}

l_b = int(raw_input("Enter l_b: "))

c = []

for i in range(0,max(a)+1):
    for j in range(0,max(a)+1):
        l = nx.shortest_path_length(G, i, j)
        if (l >= l_b):
            c.append((i, j))

G_.add_edges_from(c)

def color_or_not(node, color):
    for neighbor in G_.neighbors(node):
        neighbor_color = nodes_colors.get(neighbor, None)
        if neighbor_color == color:
            return False
    return True

def node_color(node):
    for color in colors:
        if color_or_not(node, color):
            return color

def main():
    i = 0
    for node in G_.nodes():
        nodes_colors[node] = node_color(node)
        if (nodes_colors[node] > i):
            i = nodes_colors[node]

    print math.log(i)/math.log(l_b)

main()
