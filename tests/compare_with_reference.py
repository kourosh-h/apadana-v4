"""
Compare the files of a new build with a reference build, file by file.

Usage: python compare_with_reference.py <reference_folder> <new_folder>

The header of the foyer .top file holds the date and the command line of the run;
those comment lines are skipped. All other lines must be identical.
"""
import sys
from pathlib import Path


def content_lines(path):
    lines = path.read_text().splitlines()
    if path.suffix == '.top':
        lines = [l for l in lines if not (l.startswith(';   At date') or l.startswith(';     ')
                                          or l.startswith(';   Created by') or l.startswith(';   ParmEd:') or l.startswith(';   Executable:'))]
    return lines


def main(reference, new):
    reference, new = Path(reference), Path(new)
    all_same = True
    for ref_file in sorted(p for p in reference.iterdir() if p.is_file()):
        new_file = new / ref_file.name
        if not new_file.exists():
            print(f'missing  {ref_file.name}')
            all_same = False
        elif content_lines(ref_file) == content_lines(new_file):
            print(f'same     {ref_file.name}')
        else:
            print(f'DIFFERS  {ref_file.name}')
            all_same = False
    extra = sorted(p.name for p in new.iterdir() if not (reference / p.name).exists())
    if extra:
        print('new only:', ', '.join(extra))
    return all_same


if __name__ == '__main__':
    sys.exit(0 if main(sys.argv[1], sys.argv[2]) else 1)
