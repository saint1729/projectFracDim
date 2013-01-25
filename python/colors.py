import networkx as nx
import matplotlib.pyplot as plt
G=nx.Graph()

colors = ['Red', 'Blue', 'Green', 'Yellow',  'Black', 'Pink', 'Orange', 'White', 'Gray', 'Purple', 'Brown', 'Navy']

G.add_nodes_from([1,2,3,4,5])
G.add_edges_from([(1,5),(1,3),(1,2),(1,4),(4,5)])
colors_of_nodes={}


def coloring(node, color):
   for neighbor in G.neighbors(node):
       color_of_neighbor = colors_of_nodes.get(neighbor, None)
       if color_of_neighbor == color:
          return False

   return True

def get_color_for_node(node):
    for color in colors:
       if coloring(node, color):
          return color

def main():
    for node in G.nodes():
        colors_of_nodes[node] = get_color_for_node(node)

    print colors_of_nodes


main()
