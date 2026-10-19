from setuptools import setup, find_packages

setup(
    name='apadana',
    version='4.0',
    description='Builds alkylsilane SAMs on a silica slab and writes the GROMACS input files',
    author='Kourosh Hasheminejad',
    license='MIT',
    packages=find_packages(exclude=['tests', 'examples']),
    package_data={'apadana': ['data/silica/*', 'data/forcefield/*.xml', 'data/forcefield/ForceField/*/*', 'data/mdp/*']},
    include_package_data=True,
    install_requires=[
        'mbuild==0.16.0',
        'foyer==0.12.1',
        'setuptools==69.5.1',
        'parmed',
        'matplotlib',
        'numpy',
    ],
    python_requires='>=3.9, <3.11',
    # the command-line tool: apadana build ..., apadana sheet ...
    entry_points={'console_scripts': ['apadana = apadana.command_line:main']},
)
