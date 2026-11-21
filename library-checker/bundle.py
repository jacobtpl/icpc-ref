#!/usr/bin/env python3
"""Print a solution with its local #include "..." lines expanded, ready to submit."""
import os, re, sys

seen = set()
def expand(path):
    path = os.path.realpath(path)
    if path in seen: return
    seen.add(path)
    for line in open(path):
        m = re.match(r'\s*#include "(.*)"', line)
        if m: expand(os.path.join(os.path.dirname(path), m.group(1)))
        elif line.strip() != '#pragma once': sys.stdout.write(line)

expand(sys.argv[1])
