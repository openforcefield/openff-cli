OpenFF CLI
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/openforcefield/openff-cli/workflows/ci/badge.svg)](https://github.com/openforcefield/openff-cli/actions?query=branch%3Amaster+workflow%3Aci)
[![codecov](https://codecov.io/gh/openforcefield/openff-cli/branch/master/graph/badge.svg)](https://codecov.io/gh/openforcefield/openff-cli/branch/master)


Command line utilities for the Open Force Field software stack

**Please note that this software in an early and experimental state and unsuitable for production.**

# Installation

Currently, the only available installation method is a local build

```shell
git clone https://github.com/openforcefield/openff-cli
cd openff-cli
python -m pip install .
```

# Usage

## `check_versions.py`

A simple utility to print the OpenFF packages found to be installed, and their versions

```shell
$ python openff/cli/check_versions.py
Found the following packages installed
Package name        	Version
------------            -------
OpenFF Toolkit      	0.7.2
Openforcefields     	1.2.1
OpenFF Evaluator    	Not found
OpenFF System       	Not found
OpenFF CLI          	Not found
CMILES              	0.1.5
```

## `generate_conformers.py`

Generate conformers from a starting structure, and minimize them using an OpenFF force field. Toolkit wrappers in the OpenFF Toolkit is used to call either [The RDKit](https://open-forcefield-toolkit.readthedocs.io/en/0.7.2/api/generated/openforcefield.utils.toolkits.RDKitToolkitWrapper.html#openforcefield.utils.toolkits.RDKitToolkitWrapper) or [OpenEye Omega](https://open-forcefield-toolkit.readthedocs.io/en/0.7.2/api/generated/openforcefield.utils.toolkits.OpenEyeToolkitWrapper.html#openforcefield.utils.toolkits.OpenEyeToolkitWrapper) to generate conformers, which are then energy-minimized with OpenMM.

```shell
$ python openff/cli/generate_conformers.py --help
usage: generate_conformers.py [-h] -t TOOLKIT -f FORCEFIELD -m MOLECULE
                              [-r RMS_CUTOFF] [-p PREFIX]
                              [--constrained CONSTRAINED]

Generate conformers with cheminformatics toolkits

required arguments:
  -t TOOLKIT, --toolkit TOOLKIT
                        Name of the underlying cheminformatics toolkit to use.
                        Accepted values are openeye and rdkit
  -f FORCEFIELD, --forcefield FORCEFIELD
                        Name of the force field to use, i.e. openff-1.0.0
  -m MOLECULE, --molecule MOLECULE
                        Path to an input file containing a molecule(s)

optional arguments:
  -r RMS_CUTOFF, --rms-cutoff RMS_CUTOFF
                        The redundancy cutoff between pre-minimized conformers
  -p PREFIX, --prefix PREFIX
                        The prefix for filenames of output molecules
  --constrained CONSTRAINED
                        Whether or not to use a constrained version of the
                        force field
```

### Examples:

Generate conformers from a molecule in an SDF file, using The RDKit to sample conformers, and minimizing with OpenFF 1.2.0 "Parsley":

```shell
$ python openff/cli/generate_conformers.py --molecule x.sdf --toolkit rdkit --forcefield openff-1.2.0
```

This produces many files named `molecule_N.sdf` where N begins at 0 for the lowest-energy conformer and increases with increasing conformer energy. In this case, 97 conformers were generated:

```
$ ls molecule_*.sdf | wc -l
```

Do the same, but with OpenEye Omega, OpenFF 1.0 "Parsley," and a save the conformers in files starting with `coolmol`:

```shell
$ python openff/cli/generate_conformers.py --molecule molecule.sdf --toolkit openeye --forcefield openff-1.0.0 --prefix coolmol
```

In this case, 101 conformers were generated, following the same naming scheme:

```shell
$ ls coolmol_*.sdf | wc -l
     101
```

### Copyright

Copyright (c) 2020, Open Force Field Initiative


#### Acknowledgements

Project based on the [Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.3.
