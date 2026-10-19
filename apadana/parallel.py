"""
Building several systems at the same time, one system per CPU core.

A single build runs on one core: most of its time is spent by foyer, which assigns the atom types of
the whole system in one pass. A series of systems (for example several chain lengths and coverages)
is therefore built in parallel by building every system in a process of its own. Each process runs
build_sam_system exactly as a separate run would, so the files are the same as when the systems are
built one after the other.

The processes are started with the 'spawn' method on every operating system, so that every build
starts in a fresh Python process with its own random state. A script that calls build_many_systems
therefore starts the builds under   if __name__ == '__main__':   (see examples/build_many_systems.py).
"""
import contextlib
import multiprocessing
import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

from .builder import build_sam_system

# File in the folder of every system that receives the messages of its build.
BUILD_LOG = 'build.log'


def build_one_system(settings):
    """
    Build one system with build_sam_system and write all its messages to build.log in the folder of
    the system. A failed build writes its error to the same file. Used by build_many_systems, which
    runs every call in a process of its own.

    Returns (output folder, True if the build finished, build time in minutes).
    """
    output_folder = settings['output_folder']
    os.makedirs(output_folder, exist_ok=True)
    start = time.time()
    finished = False
    # buffering=1 writes every line at once, so that a running build can be followed in its log
    with open(os.path.join(output_folder, BUILD_LOG), 'w', buffering=1) as log:
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            try:
                build_sam_system(**settings)
                finished = True
            except Exception:
                traceback.print_exc()
    minutes = (time.time() - start) / 60
    return output_folder, finished, minutes


def build_many_systems(systems, processes=None, **common_settings):
    """
    Build several systems at the same time, one system per process.

    systems         : list with one dictionary of build_sam_system arguments per system, for example
                      [{'chain_length': 6, 'coverage': 1.65}, {'chain_length': 8, 'coverage': 2.40}].
                      A system without 'output_folder' is written to L<chain_length>_G<coverage>.
    processes       : number of systems built at the same time; None gives one process per CPU core
                      (and never more processes than systems)
    common_settings : build_sam_system arguments shared by all systems, for example random_seed=1 or
                      hpc={'account': 'my_account'}; a value given in the dictionary of a system is
                      used for that system instead

    Returns a list in the order of the systems, with (output folder, True if the build finished,
    build time in minutes) for each. The messages of every build are in build.log in its folder.
    """
    # the arguments of every system: the common settings, then those given for the system
    all_settings = []
    for system in systems:
        settings = dict(common_settings)
        settings.update(system)
        if 'output_folder' not in settings:
            settings['output_folder'] = f"L{settings['chain_length']}_G{settings['coverage']}"
        # full paths, so that every process writes to the intended folder
        settings['output_folder'] = os.path.abspath(settings['output_folder'])
        if settings.get('user_surface') is not None:
            settings['user_surface'] = os.path.abspath(settings['user_surface'])
        all_settings.append(settings)

    folders = []
    for settings in all_settings:
        folders.append(settings['output_folder'])
    if len(set(folders)) < len(folders):
        raise ValueError('two systems have the same output folder; give every system its own output_folder')

    if processes is None:
        processes = os.cpu_count() or 1
    processes = max(1, min(processes, len(all_settings)))
    # flush=True: the progress is shown at once, also when the output goes to a file (for example in a batch job)
    print(f'Building {len(all_settings)} systems, {processes} at the same time. '
          f'The messages of every build are written to {BUILD_LOG} in its folder.', flush=True)

    results_by_folder = {}
    context = multiprocessing.get_context('spawn')
    with ProcessPoolExecutor(max_workers=processes, mp_context=context) as pool:
        # every submitted build, with the folder it writes to
        folder_of_build = {}
        for settings in all_settings:
            folder_of_build[pool.submit(build_one_system, settings)] = settings['output_folder']

        for n, build in enumerate(as_completed(folder_of_build), start=1):
            output_folder = folder_of_build[build]
            try:
                output_folder, finished, minutes = build.result()
            except Exception as error:
                # the process of this build stopped (for example when the computer ran out of memory)
                finished, minutes = False, float('nan')
                print(f'    process of {os.path.basename(output_folder)} stopped: {error}', flush=True)
            if finished:
                status = 'done'
            else:
                status = f'FAILED, see {BUILD_LOG}'
            print(f'[{n}/{len(all_settings)}] {os.path.basename(output_folder)}: {status} ({minutes:.1f} min)',
                  flush=True)
            results_by_folder[output_folder] = (output_folder, finished, minutes)

    results = []
    for folder in folders:
        results.append(results_by_folder[folder])
    n_finished = 0
    for output_folder, finished, minutes in results:
        if finished:
            n_finished += 1
    print(f'{n_finished} of {len(results)} systems built.', flush=True)
    return results
