#!/usr/bin/env python3

from pyprojroot import here
from pathlib import Path
import tskit, pyslim
import json
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def get_TwoSpecies_discreteWF_dirs(N_host, N_pathogen, c, s, location="slim/trees/"):
    sims = Path(here(location))
    pattern = f"TwoSpecies_discrete_WF_Nhost{N_host}_N0pathogen{N_pathogen}_c{c}_s{s}_*"
    #print(f"Pattern: {pattern}")
    dirs = sims.glob(pattern)
    return(dirs)

def read_host_pathogen_treeSeqs(dirs):
    simulations = []
    for d in dirs:
        #print(d)
        hostFile = next(d.glob("HostTreeSeq_*.trees"))
        pathogenFile = next(d.glob("PathogenTreeSeq_*.trees"))
        #print(f"hostFile: {hostFile}")
        #print(f"pathogenFile: {pathogenFile}")

        host_treeSeq = tskit.load(hostFile)
        pathogen_treeSeq = tskit.load(pathogenFile)

        if pyslim.has_vacant_samples(host_treeSeq):
            # print("Host vacant samples removed...")
            host_treeSeq = pyslim.remove_vacant(host_treeSeq).simplify()

        if pyslim.has_vacant_samples(pathogen_treeSeq):
            # print("Pathogen vacant samples removed...")
            pathogen_treeSeq = pyslim.remove_vacant(pathogen_treeSeq).simplify()


        simulations.append((host_treeSeq, pathogen_treeSeq))
    return(simulations)

############################## 
N_host = 1000
N_pathogen = 1
cValues = [0.01, 0.1, 0.5]
sValues = [0.0, 0.01, 0.1, 0.5]
# sValues = [0.01, 0.1]


rows = []

for c in cValues:
    for s in sValues: 
        print(f"Simulations: c = {c}, s = {s}")
        dirs = get_TwoSpecies_discreteWF_dirs(N_host, N_pathogen, c, s, location="slim/trees/discrete_WF")

        simulations = read_host_pathogen_treeSeqs(dirs)
        print(f"Number of simulation reps: {len(simulations)}")


        for rep, (host_ts, pathogen_ts) in enumerate(simulations):
            host_tmrcas = [tree.time(tree.root) for tree in host_ts.trees()][0]
            pathogen_tmrcas = [tree.time(tree.root) for tree in pathogen_ts.trees()][0]

            rows.append({"c": c, "s": s, "species": "host", "tmrca": host_tmrcas})
            rows.append({"c": c, "s": s, "species": "pathogen", "tmrca": pathogen_tmrcas})

df = pd.DataFrame(rows)
print(df)

#############################
# Data structure
df_wide = df.pivot_table(
    index=["c", "s"],
    columns="species",
    values="tmrca"
).reset_index()

print(df_wide)
#############################
# Plotting
 
# plt.figure()
# # sns.boxplot(data=df, x="c", y="tmrca", hue="species")
# g= sns.catplot(data=df,
#             x="c", y="tmrca", hue="species", col="s",
#             kind="box", sharey=True, showfliers=False)
# for ax in g.axes.flat:
#     ax.grid(True)

# plt.savefig(f"figures/discreteWF_HostvPathogenTMRCA_c{cValues}_s{sValues}_reps{len(simulations)}.png", dpi=300)
# plt.show()
 
plt.figure()

# g = sns.FacetGrid(df_wide, row="c", col="s", margin_titles=True)
# g.map_dataframe(sns.scatterplot, x="host", y="pathogen")

# # add grids
# for ax in g.axes.flat:
#     ax.grid(True, linestyle="--", alpha=0.5)

sns.scatterplot(data=df_wide[df_wide["s"] == 0.0], x="host", y="pathogen", hue="c")
plt.show()
