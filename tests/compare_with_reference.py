"""
Compare the files of a new build with a reference build, file by file, in all subfolders.

Usage: python compare_with_reference.py <reference_folder> <new_folder>

The comment lines at the top of a .top file (date, program and command line of the build) are
skipped. The overview sheet (.pdf, .png) is skipped too, because its files carry the date
of the build. All other files of the reference must be identical in the new build.
"""
import sys
from pathlib import Path

# Files that change with every build, also when the system is the same.
SKIPPED_SUFFIXES = ['.pdf', '.png']


def content_lines(path):
    """Lines of a file; for a .top file without the comment lines at its top."""
    lines = path.read_text().splitlines()
    if path.suffix == '.top':
        first = 0
        while first < len(lines) and lines[first].startswith(';'):
            first += 1
        lines = lines[first:]
    return lines


def compare_folders(reference, new):
    """Compare every file of the reference folder with the same file of the new folder. Returns True if all agree."""
    reference, new = Path(reference), Path(new)
    all_same = True
    n_same = 0
    n_skipped = 0
    for ref_file in sorted(p for p in reference.rglob('*') if p.is_file()):
        name = ref_file.relative_to(reference)
        new_file = new / name
        if ref_file.suffix in SKIPPED_SUFFIXES:
            n_skipped += 1
        elif not new_file.exists():
            print(f'missing  {name}')
            all_same = False
        elif content_lines(ref_file) == content_lines(new_file):
            n_same += 1
        else:
            print(f'DIFFERS  {name}')
            all_same = False

    extra = sorted(str(p.relative_to(new)) for p in new.rglob('*') if p.is_file()
                   and not (reference / p.relative_to(new)).exists())
    print(f'{n_same} files identical, {n_skipped} skipped (overview sheet)')
    if extra:
        print(f'{len(extra)} files only in the new build:', ', '.join(extra[:10]))
    return all_same


if __name__ == '__main__':
    sys.exit(0 if compare_folders(sys.argv[1], sys.argv[2]) else 1)
