import argparse
from copy import deepcopy
from typing import List, Optional

import numpy as np
from openforcefield.topology.molecule import Molecule, UndefinedStereochemistryError
from openforcefield.utils.toolkits import ToolkitRegistry
from simtk import unit

from openff.cli.core import (
    _build_simulation,
    _get_conformer_data,
    _get_forcefield,
    _minimize_conformer,
    make_registry,
)
from openff.cli.utils.utils import _enforce_dependency_version


def generate_conformers(
    molecule: str,
    registry: ToolkitRegistry,
    forcefield: str,
    constrained: bool = False,
    prefix: Optional[str] = None,
) -> List[Molecule]:

    _enforce_dependency_version("openforcefield", "0.7.1.")

    ff = _get_forcefield(forcefield, constrained)

    file_format = molecule.split(".")[-1]

    # TODO: This may not preserve order of loading molecules in
    ambiguous_stereochemistry = False
    try:
        raw_mols = registry.call(
            "from_file",
            molecule,
            file_format=file_format,
        )
    except UndefinedStereochemistryError:
        ambiguous_stereochemistry = True
        raw_mols = registry.call(
            "from_file",
            molecule,
            file_format=file_format,
            allow_undefined_stereo=True,
        )

    # When failing to parse molecules (i.e. attempting to read MOL2 with
    # RDKit, which is not supported) the toolkit can return an empty
    # list instead of raising a specific exception
    if raw_mols == []:
        from openff.cli.utils.exceptions import MoleculeParsingError

        raise MoleculeParsingError(toolkit_registry=registry, filename=molecule)

    mols = []
    for i, mol in enumerate(raw_mols):
        if prefix is not None:
            mol.name = prefix
        elif not mol.name:
            mol.name = "molecule"
        if len(raw_mols) > 1:
            mol.name += str(i)
        mols.append(mol)

    mols = _collapse_conformers(mols)

    # TODO: How to handle names of different stereoisomers? Just act like they're different conformers?
    if ambiguous_stereochemistry:
        mols_with_unpacked_stereoisomers = []
        for mol in mols:
            # TODO: This is a brute-force approach, it would be better to check stereo
            #  without needing to call enumerate_stereoisomers
            stereoisomers = mol.enumerate_stereoisomers()
            if stereoisomers:
                for i, iso in enumerate(stereoisomers):
                    iso.name = mol.name + "_stereoisomer" + str(i)
                    mols_with_unpacked_stereoisomers.append(iso)
            else:
                mols_with_unpacked_stereoisomers.append(mol)
        mols = mols_with_unpacked_stereoisomers

    for mol in mols:
        existing_conf = None
        if mol.conformers is not None:
            existing_conf = deepcopy(mol.conformers[0])
        mol.generate_conformers(
            toolkit_registry=registry,
            n_conformers=100,
            rms_cutoff=0.25 * unit.angstrom,
        )
        if existing_conf is not None:
            mol.add_conformer(existing_conf)

    # TODO: What happens if some molecules in a multi-molecule file have charges, others don't?
    mols_with_charges = []
    for mol in mols:
        if mol.partial_charges is not None:
            mols_with_charges.append(mol)

    mols_out = []
    for mol in mols:
        if mol in mols_with_charges:
            mol_with_charge = [mol]
        else:
            mol_with_charge = []
        simulation, partial_charges = _build_simulation(
            molecule=mol,
            forcefield=ff,
            mols_with_charge=mol_with_charge,
        )
        mol._partial_charges = partial_charges

        for i, conformer in enumerate(mol.conformers):
            simulation = _minimize_conformer(simulation, conformer)
            energy, positions = _get_conformer_data(simulation)
            mol = _reconstruct_mol_from_conformer(mol, positions)
            _add_metadata_to_mol(mol, energy, registry, forcefield)
            mols_out.append(mol)

    mols_out = _sort_mols(mols_out)

    return mols_out


def _collapse_conformers(molecules):
    """
    Collapse conformers of isomorphic molecules into single Molecule objects.

    This is useful when reading in multi-conformer SDF files because the SDF
    reader does not automatically collapse conformers. This function should
    not modify lists of Molecule objects whose conformers are already collapsed.

    Parameters
    ----------
    molecules : list of Molecule
        List of Molecule objects

    Returns
    -------
    collapsed_molecules : list of Molecule
        List of Molecule objects with only one object per isomorphic molecule
        and conformers of isomorphic molecules stored in each

    """
    collapsed_molecules = [molecules[0]]
    for molecule in molecules[1:]:
        if molecule == collapsed_molecules[-1]:
            for conformer in molecule.conformers:
                collapsed_molecules[-1].add_conformer(conformer)
        else:
            collapsed_molecules.append(molecule)

    return collapsed_molecules


def _reconstruct_mol_from_conformer(
    mol: Molecule,
    positions: unit.Quantity,
) -> Molecule:
    mol = deepcopy(mol)
    mol._conformers = None
    min_coords = (
        np.array([[atom.x, atom.y, atom.z] for atom in positions]) * unit.nanometer
    )
    mol.add_conformer(min_coords)
    return mol


def _add_metadata_to_mol(
    mol: Molecule,
    energy: unit.Quantity,
    toolkit_registry: ToolkitRegistry,
    ff_name: str,
) -> None:
    toolkit_name = [*toolkit_registry.registered_toolkit_versions][0]
    toolkit_version = toolkit_registry.registered_toolkit_versions[toolkit_name]
    mol.properties["absolute energy (kcal/mol): "] = energy
    mol.properties["conformer generation toolkit: "] = (
        toolkit_name + " " + toolkit_version
    )
    mol.properties["minimized against: "] = ff_name


def _sort_mols(mols: List[Molecule]) -> List[Molecule]:
    final_list = []
    unique_mol_names = set([mol.name for mol in mols])

    # To handle multi-molecule SDF files, group molecules by molecules, and then conformers
    for mol_name in unique_mol_names:
        mols_with_energy = [
            (mol, mol.properties["absolute energy (kcal/mol): "])
            for mol in mols
            if mol.name == mol_name
        ]
        mols_with_energy.sort(key=lambda x: x[1])
        # TODO: Here could be some logic for filtering out high-energy conformers
        sorted_mols = [x[0] for x in mols_with_energy]
        for i, sorted_mol in enumerate(sorted_mols):
            sorted_mol.name += "_conf" + str(i)
            final_list.append(sorted_mol)
    return final_list


def write_mols(
    mols: List[Molecule],
    toolkit_registry: ToolkitRegistry,
    molecule: str,
    prefix: Optional[str],
) -> None:
    """Save minimized structures, with data in SD tags, to files"""
    if prefix is not None:
        prefix_out = prefix
    else:
        prefix_from_filename = ".".join(molecule.split(".")[:-1])
        prefix_out = prefix_from_filename

    if len(mols) == 1:
        mols[0].to_file(
            file_path=prefix_out + "_0.sdf",
            file_format="SDF",
            toolkit_registry=toolkit_registry,
        )
    else:
        for i, mol in enumerate(mols):
            filename = prefix_out + "_" + str(i) + ".sdf"
            mol.to_file(
                file_path=filename, file_format="SDF", toolkit_registry=toolkit_registry
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate conformers with cheminformatics toolkits"
    )
    parser._action_groups.pop()
    required_args = parser.add_argument_group("required arguments")
    optional_args = parser.add_argument_group("optional arguments")

    required_args.add_argument(
        "-t",
        "--toolkit",
        type=str,
        required=True,
        help="Name of the underlying cheminformatics toolkit to use. Accepted"
        " values are openeye and rdkit",
    )
    required_args.add_argument(
        "-f",
        "--forcefield",
        type=str,
        required=True,
        help="Name of the force field to use, i.e. openff-1.0.0",
    )
    required_args.add_argument(
        "-m",
        "--molecule",
        type=str,
        required=True,
        help="Path to an input file containing a molecule(s)",
    )
    optional_args.add_argument(
        "-r",
        "--rms-cutoff",
        type=float,
        default=0.25,
        help="The redundancy cutoff between pre-minimized conformers",
    )
    optional_args.add_argument(
        "-p",
        "--prefix",
        type=str,
        help="The prefix for filenames of output molecules",
    )
    optional_args.add_argument(
        "--constrained",
        type=bool,
        default=False,
        help="Whether or not to use a constrained version of the force field",
    )
    args = parser.parse_args()

    registry = make_registry(args.toolkit)

    mols = generate_conformers(
        molecule=args.molecule,
        registry=registry,
        forcefield=args.forcefield,
        constrained=args.constrained,
        prefix=args.prefix,
    )

    write_mols(
        mols, toolkit_registry=registry, molecule=args.molecule, prefix=args.prefix
    )
