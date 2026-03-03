#!/usr/bin/env python3

import msprime
import tskit
from IPython.display import display


# Parameters (Host)
num_samples = 5
Ne = 10000  # Effective population size
L = 1e5
r = 1e-9

### Host simulation
## Simulate a coalescent tree sequence with recombination
ts = msprime.sim_ancestry(
    samples=num_samples,
    population_size    = Ne,
    recombination_rate = r,
    sequence_length    = L,
    ploidy             = 1)

ts.dump("host_treeSequence.trees")

ts_host = tskit.load("host_treeSequence.trees")


# Host tree sequence
css = """
.node:not(.leaf) > .sym, .node:not(.leaf) > .lab {display: none}
#myUID .background * {fill: #6666FF}
.edge {stroke: darkblue; stroke-width: 4px;}
"""


svg_size   = (1200, 250)
y_tick_pos = [0, 1000, 5000, 10000]

svg_ts_host = ts_host.draw_svg(size = svg_size, canvas_size = (1400, 250),
                               y_axis=True,
                               #x_scale = "treewise",
                               root_svg_attributes = {'id': "myUID"},
                               style = css)


with open('host_treeSequence.svg', 'w') as f:
    f.write(svg_ts_host)

