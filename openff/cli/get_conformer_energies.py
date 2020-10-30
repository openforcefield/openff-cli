import argparse
from typing import List

from openforcefield.topology import Molecule
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


def get_conformer_energies(
    molecule: str,
    registry: ToolkitRegistry,
    forcefield: str,
    constrained: bool = False,
) -> List[Molecule]:

    _enforce_dependency_version("openforcefield", "0.7.0")

    file_format = molecule.split(".")[-1]

    loaded_molecules = registry.call(
        "from_file",
        molecule,
        file_format=file_format,
    )

    if type(loaded_molecules) is not list:
        loaded_molecules = [loaded_molecules]

    mols = [loaded_molecules[0]]
    for mol in loaded_molecules[1:]:
        if mol == mols[-1]:
            for conformer in mol.conformers:
                mols[-1].add_conformer(conformer)
        else:
            mols.append(molecule)

    n_molecules = len(mols)
    n_conformers = sum([mol.n_conformers for mol in mols])
    print(
        f"{n_molecules} unique molecule(s) loaded, with {n_conformers} total conformers"
    )

    ff = _get_forcefield(forcefield, constrained)

    mols_with_charges = []
    for mol in mols:
        if mol.partial_charges is not None:
            mols_with_charges.append(mol)

    # This is duplicated from generate_conformers
    minimized_mols = []
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

        mol.properties["minimized against: "] = forcefield
        mol.properties["original conformer energies (kcal/mol)"] = mol.n_conformers * [
            None
        ]
        mol.properties["minimized conformer energies (kcal/mol)"] = mol.n_conformers * [
            None
        ]
        for i, conformer in enumerate(mol.conformers):
            simulation.context.setPositions(conformer)
            pre_energy, pre_positions = _get_conformer_data(simulation)
            mol.properties["original conformer energies (kcal/mol)"][i] = pre_energy

            simulation = _minimize_conformer(simulation, conformer)
            min_energy, min_positions = _get_conformer_data(simulation)
            mol.properties["minimized conformer energies (kcal/mol)"][i] = min_energy
            mol.conformers[i] = min_positions
        minimized_mols.append(mol)

    return minimized_mols


def _print_mol_data(mols):
    pre_key = "original conformer energies (kcal/mol)"
    min_key = "minimized conformer energies (kcal/mol)"
    for mol_idx, mol in enumerate(mols):
        forcefield = mol.properties["minimized against: "]
        print(f"printing mol {mol_idx}, minimized against {forcefield}")
        for conformer_idx in range(mol.n_conformers):
            pre_energy = mol.properties[pre_key][conformer_idx]
            min_energy = mol.properties[min_key][conformer_idx]
            print(
                "%5d / %5d : %8.3f kcal/mol %8.3f kcal/mol  %8.3f Angstroms"
                % (
                    conformer_idx + 1,
                    mol.n_conformers,
                    pre_energy / unit.kilocalories_per_mole,
                    min_energy / unit.kilocalories_per_mole,
                    0.0,  # RMSD goes here
                )
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Evaluate conformer energies with OpenMM"
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
        help="Path to an input file containing a molecule(s), single or multi-conformers",
    )
    optional_args.add_argument(
        "--constrained",
        type=bool,
        default=False,
        help="Whether or not to use a constrained version of the force field",
    )
    args = parser.parse_args()

    registry, _ = make_registry(args.toolkit)

    mols = get_conformer_energies(
        molecule=args.molecule,
        registry=registry,
        forcefield=args.forcefield,
        constrained=args.constrained,
    )

    _print_mol_data(mols)
