#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


N, beta, G = 1000, 0.1, 20
cVs = [0.3, 0.6, 0.8, 1]
data_dir = Path(".")       



def read_txt(path):
    data = {}
    sp = None
    for line in Path(fname).read_text().splitlines():
        if line.startswith("# Size history"):
            sp = line.split("species ")[1].rstrip("):")
        elif sp and line[:1].isdigit():
            data[sp] = [int(x) for x in line.split(",")]
            sp = None

    df = pd.DataFrame(data)
    df["x"] = df["pathogen"] / df["host"]            # fraction infected
    #end = df.index[df.index % (G+1) == G]
    df["cycle"] = df.index // (G+1)          # host generation
    df["g"]     = df.index %  (G+1)          # tick within generation (0..G)
    df["t"]     = df["cycle"] + df["g"]/(G+1)  # continuous time in host cycles
    df.index.name = "tick"
    return df

### Run for params ###

def x_star(cV, beta, G):
    """Stationary infected fraction; 0 if the pathogen cannot persist."""
    e = np.exp(-beta * G)
    return max((cV - e) / (cV * (1 - e)), 0.0)
 
 
# ------------------------------------------------------------------ plotting
fig, ax = plt.subplots(figsize=(10, 4), layout="constrained")
# colors = plt.cm.viridis(np.linspace(0, 0.9, len(cVs)))
palette = "Reds"
colors = sns.color_palette(palette, n_colors=len(cVs) + 1)[1:]   # light = low cV, dark = high cV
HostCyclesN = 0.1
 
for cV, col in zip(cVs, colors):
    fname = data_dir / f"N{N}_beta{beta}_cV{cV}_G{G}.txt"
    if not fname.exists():
        print(f"missing: {fname.resolve()}")
        continue
 
    df = read_txt(fname).dropna()
    end = df[df["g"] == G]                 # x_n: end of each host generation
    print(end)
 
    xs = x_star(cV, beta, G)
    ax.plot(df["t"], df["x"], "-", lw=1.5, c=col, alpha=0.35)
    ax.plot(end["t"], end['x'], lw=2, c=col, label=rf"$c_v={cV}$")
    ax.axhline(xs, ls="--", c=col, lw=2)
    #last_cycle = max(last_cycle, end["cycle"].max() + 1)
    ax.text(1.01, xs, rf"$x^*={xs:.2f}$", c=col, va="center", fontsize=9,
            transform=ax.get_yaxis_transform())
ax.set(xlabel=r"Host generation $n$", ylabel=r"Infected fraction $x_n$",
       xlim=(0, HostCyclesN * N),
       ylim=(0, 1.05))
       #title=rf"$N={N},\ \beta={beta},\ G={G}$")

# ax.text(1.02, 0.10, fr"$N = {N}$", va="center", fontsize=10,
#             transform=ax.get_yaxis_transform())
# ax.text(1.02, 0.05, fr"$\beta = {beta}$", va="center", fontsize=10,
#             transform=ax.get_yaxis_transform())
# ax.text(1.02, 0.0, fr"$G = {G}$", va="center", fontsize=10,
#             transform=ax.get_yaxis_transform())
ax.text(0.8, 1.08, fr"$N = {N}$; $\beta = {beta}$; $G = {G}$", va='center', fontsize=10,
        transform=ax.get_yaxis_transform())

handles, labels = ax.get_legend_handles_labels()
fig.legend(handles[::-1], labels[::-1],
          loc="outside right", bbox_to_anchor=(1.03, 0.23), ncol=1, frameon=True, fontsize=8)
#ax.set_title(fr"$N = {N}$, $\beta = {beta}$, $G = {G}$", loc='left', fontsize=10)
#plt.tight_layout()
plt.grid(True)
plt.savefig(f"Demographies_N{N}_beta{beta}_G{G}_cycles{HostCyclesN}.pdf", bbox_inches='tight', dpi=300)
plt.show()

### 
# fname = f"N{N}_beta{beta}_cV{cV}_G{G}.txt"   # N100_beta0.1_cV0.8_G10.txt
# print(Path(fname))

# x_star = (cV - np.exp(-beta*G)) / (cV * (1 - np.exp(-beta*G)))

# fig, ax = plt.subplots(figsize=(8, 4))
# ax.plot(df["t"], df["x"], lw=1, label=r"$X_n$")
# ax.plot(end["t"], end["x"], "x", ms=3)
# # ax.plot(df['cycle'], df["x"], lw=1, label="infected fraction")
# # ax.plot(df['cycle'], df.loc[end, "x"], "o", ms=4, label=r"$x_n$ (end of generation)")
# ax.axhline(x_star, ls="--", c="k", label=rf"$x^*={x_star:.2f}$")
# ax.set(xlabel="tick", ylabel=r"Infected fraction $X_n$", title=f"N={N}, beta={beta}, cV={cV}, G={G}")
# ax.legend()
# plt.tight_layout()
# plt.grid(True)
# plt.show()
