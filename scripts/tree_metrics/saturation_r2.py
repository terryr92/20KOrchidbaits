#!/usr/bin/env python3
"""
Compute saturation as R^2 of the regression:
p-distance (from alignment) vs patristic distance (from tree)

Inputs:
- phylip alignment (same used for IQ-TREE)
- newick treefile

Encoding:
- '0' and '1' are states
- '?' is missing and ignored pairwise
"""
import sys
import numpy as np
from ete3 import Tree

phy = sys.argv[1]
treefile = sys.argv[2]

# Read PHYLIP
with open(phy, "r") as f:
    header = f.readline().split()
    n = int(header[0]); L = int(header[1])
    names, seqs = [], []
    for _ in range(n):
        line = f.readline().strip().split()
        names.append(line[0])
        seqs.append(line[1])

seqs = dict(zip(names, seqs))

t = Tree(treefile, format=1)
tips = [leaf.name for leaf in t.get_leaves()]
tips = [x for x in tips if x in seqs]

def p_dist(a, b):
    sa, sb = seqs[a], seqs[b]
    diff = 0; comp = 0
    for x, y in zip(sa, sb):
        if x == "?" or y == "?":
            continue
        comp += 1
        diff += int(x != y)
    return (diff / comp) if comp else np.nan

xs = []
ys = []
for i in range(len(tips)):
    for j in range(i+1, len(tips)):
        a, b = tips[i], tips[j]
        pd = p_dist(a, b)
        if np.isnan(pd):
            continue
        pat = t.get_distance(a, b)
        xs.append(pat)
        ys.append(pd)

x = np.array(xs); y = np.array(ys)
if len(x) < 3:
    print("nan")
    sys.exit(0)

# Linear regression y = m*x + c; compute R^2
m, c = np.polyfit(x, y, 1)
yhat = m*x + c
ss_res = np.sum((y - yhat)**2)
ss_tot = np.sum((y - np.mean(y))**2)
r2 = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
print(r2)
