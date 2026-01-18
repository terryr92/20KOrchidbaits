#!/usr/bin/env python3
"""
Extract median support from an IQ-TREE Newick file.
This assumes supports are stored as node labels (common IQ-TREE output).
"""
import re
import statistics
import sys

treefile = sys.argv[1]
txt = open(treefile, "r").read()

# Common pattern: )95: or )95.3:
supports = [float(x) for x in re.findall(r"\)(\d+(?:\.\d+)?)\s*:", txt)]
med = statistics.median(supports) if supports else float("nan")
print(med)
