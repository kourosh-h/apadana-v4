"""
Build a series of systems at the same time, one system per CPU core.

Here: the chain lengths L = 4 ... 8 at the coverage 2.40 molecules/nm^2, with random seed 1. Every
system is written to its own folder (L4_G2.4, L5_G2.4, ...), and the messages of its build to
build.log in that folder. A build needs between 0.5 and 0.8 GB of memory, so the number of builds that
run at the same time (processes) is chosen to fit the free memory of the computer.

Run:  python build_many_systems.py
"""
from apadana import build_many_systems

# The builds start in new Python processes, which read this file again;
# the line below keeps them from starting the builds a second time.
if __name__ == '__main__':
    systems = []
    for chain_length in [4, 5, 6, 7, 8]:
        systems.append({'chain_length': chain_length, 'coverage': 2.40})
    build_many_systems(systems, processes=4, random_seed=1)
