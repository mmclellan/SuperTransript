import sys
from networkx import *

print(sys.version)

print(networkx.__version__)

G = Graph()

G.add_node("spam")
G.add_edge(1,2)
print(G.nodes())

print(G.edges())

