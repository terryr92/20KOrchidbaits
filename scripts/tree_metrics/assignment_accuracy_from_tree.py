#!/usr/bin/env python3
"""
Compute assignment accuracy from a tree by checking whether each pollinarium
clusters within its expected species clade.
This requires:
- a mapping file with columns: sample,species,type
  where type is 'pollinarium' or 'reference'
Rule (simple, reproducible):
- For each pollinarium, find its closest reference (by path length) and count as correct
  if closest reference species matches expected species.
"""
import sys
import pandas as pd
from ete3 import Tree

treefile = sys.argv[1]
mapping  = sys.argv[2]
out_csv  = sys.argv[3]

meta = pd.read_csv(mapping)
t = Tree(treefile, format=1)

# Precompute references
refs = meta[meta["type"] == "reference"].copy()
poll = meta[meta["type"] == "pollinarium"].copy()

ref_species = dict(zip(refs["sample"], refs["species"]))
poll_species = dict(zip(poll["sample"], poll["species"]))

# Ensure tips exist
tips = set([leaf.name for leaf in t.get_leaves()])
refs = refs[refs["sample"].isin(tips)]
poll = poll[poll["sample"].isin(tips)]

def closest_reference(poll_tip):
    best = None
    best_d = None
    for r in refs["sample"]:
        d = t.get_distance(poll_tip, r)
        if best_d is None or d < best_d:
            best_d = d
            best = r
    return best, best_d

rows = []
correct = 0
for p in poll["sample"]:
    r, d = closest_reference(p)
    exp = poll_species[p]
    pred = ref_species[r]
    ok = (exp == pred)
    correct += int(ok)
    rows.append({"pollinarium": p, "expected": exp, "closest_ref": r, "predicted": pred, "distance": d, "correct": ok})

acc = correct / len(rows) if rows else float("nan")
pd.DataFrame(rows).to_csv(out_csv, index=False)
print(acc)
