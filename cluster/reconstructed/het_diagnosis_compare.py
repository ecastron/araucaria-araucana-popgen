#!/usr/bin/env python3
"""Compare a reconstructed het_diagnosis.json with the original, value by value.

Usage: het_diagnosis_compare.py <original.json> <reconstructed.json>
Prints every leaf that differs (strings exactly, numbers with their absolute difference) and the
max absolute numeric difference per top-level block. Written 18-19 Nov 2026 alongside
het_diagnosis_assemble.py; this is a checker, it produces no project file.
"""
import json, sys
A, B = (json.load(open(p)) for p in sys.argv[1:3])
diffs, maxd = [], {}
def walk(a, b, path, block):
    if isinstance(a, dict) or isinstance(b, dict):
        for k in sorted(set(a or {}) | set(b or {})):
            walk((a or {}).get(k), (b or {}).get(k), f"{path}.{k}", block or k)
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, f"{path}[{i}]", block)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        d = abs(a - b); maxd[block] = max(maxd.get(block, 0), d)
        if d > 0: diffs.append((path, a, b, d))
    elif a != b:
        diffs.append((path, a, b, None)); maxd.setdefault(block, 0)
    else:
        maxd.setdefault(block, 0)
walk(A, B, "", None)
for blk, d in maxd.items():
    print(f"{blk:28s} max |diff| = {d:.6g}")
print(f"\n{len(diffs)} differing leaves:")
for p, a, b, d in diffs:
    print(f"  {p}\n     original: {a!r}\n     reconstr: {b!r}" + (f"\n     |diff| = {d:.6g}" if d is not None else ""))
