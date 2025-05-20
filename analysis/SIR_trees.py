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


#print(ts)
print(ts.draw_text())


md = ts.metadata['SLiM']['user_metadata']
infection_dict = {pedigree_id: InfectionTag for pedigree_id, InfectionTag in zip(md['ID'], md['InfectionTags'])}
print(infection_dict)

susceptible = 0
resistant = 1
infected = 2


def filter_nodes(ts, SIRclass):
    subset_nodes = []
    for node in ts.nodes():
        if node.individual == -1:
            continue
        ind = ts.individual(node.individual)
        id = ind.metadata['pedigree_id']
        if id in infection_dict and infection_dict[id] == SIRclass:
            subset_nodes.append(node.id)

    print(f"Length: {len(subset_nodes)}")
    return(subset_nodes)


print("Infected:")
infected_nodes = filter_nodes(ts, infected)
print(infected_nodes)
ts_infected = ts.simplify(samples=infected_nodes)
print(ts_infected)
print(ts_infected.draw_text())



print("Susceptible: ")
susceptible_nodes = filter_nodes(ts, susceptible)
print(susceptible_nodes)
ts_susceptible = ts.simplify(samples=susceptible_nodes)
print(ts_susceptible)
print(ts_susceptible.draw_text())


print("Resistant: ")
resistant_nodes = filter_nodes(ts, resistant)
print(resistant_nodes)
ts_resistant = ts.simplify(samples=resistant_nodes)
print(ts_resistant)
print(ts_resistant.draw_text())



# Calculate TMRCA of trees

print(f"\nTree Root length (TMRCA): ")

def tmrca(ts):
    for tree in ts.trees():
        root = ts.node(tree.root)
        print(root)
        tmrca = root.time
        print(tmrca)
    return(tmrca)

tmrca_resistant = tmrca(ts_resistant)
tmrca_susceptible = tmrca(ts_susceptible)
tmrca_infected = tmrca(ts_infected)
