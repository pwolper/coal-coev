#!/usr/bin/env python3

from pyprojroot import here
from pathlib import Path
import tskit, pyslim
import json
import glob


def get_TwoSpecies_discreteWF_dirs(N_host, N_pathogen, c, location="slim/trees"):
    sims = Path(here(location))
    pattern = f"TwoSpecies_discrete_WF_Nhost{N_host}_Npathogen{N_pathogen}_c{c}_*"
    print(f"Pattern: {pattern}")

    dirs = sims.glob(pattern)
    return(dirs)

def read_host_pathogen_treeSeqs(dirs):
    simulations = []
    for d in dirs:
        print(d)
        hostFile = next(d.glob("HostTreeSeq_*.trees"))
        pathogenFile = next(d.glob("PathogenTreeSeq_*.trees"))
        #print(f"hostFile: {hostFile}")
        #print(f"pathogenFile: {pathogenFile}")

        host_treeSeq = tskit.load(hostFile)
        pathogen_treeSeq = tskit.load(pathogenFile)

        if pyslim.has_vacant_samples(host_treeSeq):
            print("Host vacant samples removed...")
            host_treeSeq = pyslim.remove_vacant(host_treeSeq)

        if pyslim.has_vacant_samples(pathogen_treeSeq):
            print("Pathogen vacant samples removed...")
            pathogen_treeSeq = pyslim.remove_vacant(pathogen_treeSeq)

        simulations.append((host_treeSeq, pathogen_treeSeq))
    return(simulations)

############################## 
N_host = 10
N_pathogen = 1
c = 0.1

dirs = get_TwoSpecies_discreteWF_dirs(N_host, N_pathogen, c)
print(dirs)

simulations = read_host_pathogen_treeSeqs(dirs)
print(simulations)


host_treeSeq, pathogen_treeSeq = simulations[0]

# for tree in host_treeSeq.trees():
#     print(tree.num_roots)

# for tree in pathogen_treeSeq.trees():
#     print(tree.num_roots)


# print(host_treeSeq.draw_text())
# print(pathogen_treeSeq.draw_text())
