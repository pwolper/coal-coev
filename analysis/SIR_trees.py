#!/usr/bin/env python3

from pyprojroot import here
import tskit
import json

ts = tskit.load(here("slim/trees/Moran_SIR_test.trees"))

valid_nodes = [node.id for node in ts.nodes() if node.metadata.get("is_vacant") == [0]]

# Simplify the tree sequence with these nodes as samples
ts = ts.simplify(samples=valid_nodes)

# for tree in ts.trees():
#     assert tree.num_roots == 1, 'Tree not coalesced!'


print(ts)
print(ts.draw_text())

md = ts.metadata['SLiM']['user_metadata']
infection_dict = {pedigree_id: InfectionTag for pedigree_id, InfectionTag in zip(md['ID'], md['InfectionTags'])}
print(infection_dict)

susceptible = 0
resistant = 1
infected = 2

subset_nodes = []

for node in ts.nodes():
    if node.individual == -1:
        continue
    ind = ts.individual(node.individual)
    id = ind.metadata['pedigree_id']
    if id in infection_dict and infection_dict[id] == infected:
        subset_nodes.append(node.id)

print(subset_nodes)

ts_infected = ts.simplify(samples=subset_nodes)

print(ts_infected)
print(ts_infected.draw_text())
