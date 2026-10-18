"""
Build the test case of apadana: chain length L = 6 and coverage 1.65 molecules/nm^2 on the relaxed
slab of the package, with random seed 1, so that every build gives the same system.

The output is written to test_case_L6_G1.65/ next to this script. The stored copy of this folder
shows every file that apadana writes for a SAM system, and tests/run_test_case.py builds it again
and compares the new files with the stored ones.

Run:  python build_test_case.py
"""
from pathlib import Path

from apadana import build_sam_system

HERE = Path(__file__).parent

build_sam_system(chain_length=6, coverage=1.65, random_seed=1, output_folder=HERE / 'test_case_L6_G1.65')
