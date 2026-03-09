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

            rows.append({"c": c, "s": s, "tmrca_host": host_tmrcas, "tmrca_pathogen": pathogen_tmrcas})

df = pd.DataFrame(rows)
print(df)

#############################
# # Data structure
df_long = df.melt(id_vars=["c", "s"],
                  value_vars=["tmrca_host", "tmrca_pathogen"],
                  var_name="species",
                  value_name="tmrca")
df_long["species"] = df_long["species"].str.replace("tmrca_", "")

print(df_long)
#############################
# Plotting
 
row_colors = {"host": "tab:blue", "pathogen": "tab:orange"}
# sns.boxplot(data=df, x="c", y="tmrca", hue="species")
g = sns.catplot(data=df_long,
               x="c", y="tmrca", hue="species", row="species", col="s",
                kind="box", sharey='row', boxprops={"alpha": 0.4}, showfliers=False, width=0.6)
# g.map_dataframe(sns.stripplot, x="c", y="tmrca", jitter=True)
# loop over rows (species) and columns (s values)
species_order = ["host", "pathogen"]
s_values = sorted(df_long["s"].unique())

for row_idx, species in enumerate(species_order):
    for col_idx, s_val in enumerate(s_values):
        ax = g.axes[row_idx, col_idx]
        # select only the subset for this facet
        df_facet = df_long[(df_long["species"] == species) & (df_long["s"] == s_val)]
        sns.stripplot(
            data=df_facet,
            x="c",
            y="tmrca",
            # jitter=0.25,
            alpha=0.8,
            color=row_colors[species],
            ax=ax, legend=False)

for ax in g.axes.flat:
    ax.grid(True)

g.fig.savefig(f"{here("analysis")}/figures/discreteWF_HostvPathogenTMRCA_c{cValues}_s{sValues}_reps{len(simulations)}.png", dpi=300)
 
# plt.figure()

# g = sns.FacetGrid(df, row="c", col="s", margin_titles=True)
# g.map_dataframe(sns.scatterplot, x="tmrca_host", y="tmrca_pathogen")

# # add grids
# for ax in g.axes.flat:
#     ax.grid(True, linestyle="--", alpha=0.5)

# # sns.scatterplot(data=df_wide[df_wide["s"] == 0.0], x="host", y="pathogen", hue="c")
# plt.show()
