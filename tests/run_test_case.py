"""
Check that apadana still builds the test case exactly as stored in examples/test_case_L6_G1.65.

The test case is built again in a temporary folder with the same settings (L = 6, coverage
1.65 molecules/nm^2, random seed 1), and every stored file is compared with the new one
(see compare_with_reference.py). The build takes about 10 minutes.

Run:  python tests/run_test_case.py
"""
import sys
import tempfile
from pathlib import Path

from apadana import build_sam_system
from compare_with_reference import compare_folders

REFERENCE = Path(__file__).parent.parent / 'examples' / 'test_case_L6_G1.65'


def main():
    with tempfile.TemporaryDirectory() as temporary_folder:
        new = Path(temporary_folder) / 'test_case_L6_G1.65'
        build_sam_system(chain_length=6, coverage=1.65, random_seed=1, output_folder=new)
        print('\nComparison with the stored test case:')
        all_same = compare_folders(REFERENCE, new)
    if all_same:
        print('TEST PASSED: the build is identical to the stored test case.')
    else:
        print('TEST FAILED: the build differs from the stored test case.')
    return all_same


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
