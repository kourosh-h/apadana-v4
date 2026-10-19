"""Run the command-line tool as  python -m apadana ...  (see command_line.py)."""
import sys

from .command_line import main

if __name__ == '__main__':
    sys.exit(main())
