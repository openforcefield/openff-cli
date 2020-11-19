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

## `openff check_versions`

A simple utility to print the OpenFF packages found to be installed, and their versions

```shell
$ openff check_versions
Found the following packages installed
Package name        	Version
------------            -------
OpenFF Toolkit      	0.8.0
Openforcefields     	1.3.0
OpenFF Evaluator    	Not found
OpenFF System       	Not found
OpenFF CLI          	0.0.0+100.ge013bfc.dirty
CMILES              	v0.1.6
```

## `openff generate_conformers`

Generate conformers from a starting structure, and minimize them using an OpenFF force field. Toolkit wrappers in the OpenFF Toolkit is used to call either [The RDKit](https://open-forcefield-toolkit.readthedocs.io/en/0.7.2/api/generated/openforcefield.utils.toolkits.RDKitToolkitWrapper.html#openforcefield.utils.toolkits.RDKitToolkitWrapper) or [OpenEye Omega](https://open-forcefield-toolkit.readthedocs.io/en/0.7.2/api/generated/openforcefield.utils.toolkits.OpenEyeToolkitWrapper.html#openforcefield.utils.toolkits.OpenEyeToolkitWrapper) to generate conformers, which are then energy-minimized with OpenMM.

This requires [`OpenFF Toolkit`](https://github.com/openforcefield/openforcefield) version 0.7.1 or newer.

```shell
$ openff generate_conformers --help
Usage: openff generate_conformers [OPTIONS]

  Generate conformers from a starting structure, and minimize them using an
  OpenFF force field.

Options:
  -t, --toolkit [openeye|rdkit]  Name of the underlying cheminformatics
                                 toolkit to use. Accepted values are openeye
                                 and rdkit  [required]

  -f, --forcefield TEXT          Name of the force field to use, i.e.
                                 openff-1.0.0  [required]

  -m, --molecule TEXT            Path to an input file containing a
                                 molecule(s)  [required]

  -r, --rms-cutoff FLOAT         The redundancy cutoff between pre-minimized
                                 conformers

  -p, --prefix TEXT              The prefix for filenames of output molecules
  --constrained BOOLEAN          Whether or not to use a constrained version
                                 of the force field

  --help                         Show this message and exit.
```

### Examples:

Generate conformers from a molecule in an SDF file, using The RDKit to sample conformers, and minimizing with OpenFF 1.2.0 "Parsley":

```shell
$ openff generate_conformers --molecule x.sdf --toolkit rdkit --forcefield openff-1.2.0
```

This produces many files named `molecule_N.sdf` where N begins at 0 for the lowest-energy conformer and increases with increasing conformer energy. In this case, 97 conformers were generated:

```
$ ls molecule_*.sdf | wc -l
```

Do the same, but with OpenEye Omega, OpenFF 1.0 "Parsley," and a save the conformers in files starting with `coolmol`:

```shell
$ openff generate_conformers --molecule molecule.sdf --toolkit openeye --forcefield openff-1.0.0 --prefix coolmol
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
