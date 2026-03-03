#!/usr/bin/env python3

import msprime
import tskit
from IPython.display import display


# Parameters (Pathogen)
num_samples = 5
Ne = 100  # Effective population size
L = 1e4
r = 1e-6

## Pathogen simulation
# Simulate a coalescent tree sequence with recombination
ts = msprime.sim_ancestry(
    samples=num_samples,
    population_size    = Ne,
    recombination_rate = r,
    sequence_length    = L,
    ploidy             = 1)

ts.dump("pathogen_treeSequence.trees")

ts_pathogen = tskit.load("pathogen_treeSequence.trees")


# Pathogen tree sequence
css = """
.node:not(.leaf) > .sym, .node:not(.leaf) > .lab {display: none}
#myUID .background * {fill: #FF6666}
.edge {stroke: darkred; stroke-width: 4px;}
"""

# .node > .sym, .node > .lab {display: none}
# .node.leaf > .sym, .node.leaf > .lab {display: initial}

svg_size   = (1200, 250)
y_tick_pos = [0, 1000, 5000, 10000]

svg_ts_pathogen = ts_pathogen.draw_svg(size=svg_size, canvas_size = (1400, 250),
                                         y_axis=True,
                                         # x_scale   = "treewise",
                                         root_svg_attributes={'id': "myUID"},
                                         style=css)


with open('pathogen_treeSequence.svg', 'w') as f:
    f.write(svg_ts_pathogen)
