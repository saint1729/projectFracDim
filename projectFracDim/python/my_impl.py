import networkx as nx

G = nx.Graph()

l_B_max = diameter(G) + 1

N = G.size()

c = [[0 for col in range(l_B_max)] for row in range(N)]

if l_B == 1 :
    N_B = N

if l_B >= l_B_max :
    N_B = 1

id = 0
for l in range(l_B_max) :
    c[id][l] = 0

id = 1

color = 1
while (i < N) :
    l = [[0 for col in range(j)] for row in range(i)]
    while (j < i) :
        l[i][j] = shortest_path_length(G, i, j)
        j = j + 1
    l_B = 1
    while (l_B <= l_B_max) :
        while (j < i) :
            if l[i][j] >= l_B :
                c[j][l[i][j]] = color + 1
            j = j + 1
        l_B = l_B + 1
    i = i + 1
