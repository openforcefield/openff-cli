import argparse
from typing import List, Optional

from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.utils.toolkits import ToolkitRegistry
from simtk import unit

from openff.cli.generate_conformers import (
    _build_simulation,
    _get_minimized_data,
    make_registry,
)
from openff.cli.utils.utils import _enforce_dependency_version


def get_conformer_energies(
    molecule: str,
    registry: ToolkitRegistry,
    forcefield: str,
    constrained: bool = False,
    prefix: Optional[str] = None,
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
    for molecule in loaded_molecules[1:]:
        if molecule == mols[-1]:
            for conformer in molecule.conformers:
                mols[-1].add_conformer(conformer)
        else:
            mols.append(molecule)

    n_molecules = len(mols)
    n_conformers = sum([mol.n_conformers for mol in mols])
    print(
        f"{n_molecules} unique molecule(s) loaded, with {n_conformers} total conformers"
    )

    # This is duplicated from generate_conformers
    ff_name = forcefield
    if constrained:
        split_name = ff_name.split("-")
        split_name[0] = "openff_unconstrained"
        ff_name = "-".join(split_name)
    if not ff_name.endswith(".offxml"):
        ff_name += ".offxml"
    ff = ForceField(ff_name)

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

        mol.properties["minimized against: "] = ff_name
        mol.properties["conformer energies (kcal/mol)"] = list()
        for i, conformer in enumerate(mol.conformers):
            energy, positions = _get_minimized_data(conformer, simulation)
            # mol = _reconstruct_mol_from_conformer(mol, positions)
            # This function could be moved to a common place and modified
            # to more generally take in a set of optional args
            mol.properties["conformer energies (kcal/mol)"].append(energy)
            mol.conformers[i] = positions
        minimized_mols.append(mol)

    _print_mol_data(mols=minimized_mols, prefix=prefix)


def _print_mol_data(mols, prefix=None):
    if prefix is None:
        prefix = "molecules"

    for mol_idx, mol in enumerate(mols):
        forcefield = mol.properties["minimized against: "]
        print(f"printing mol {mol_idx}, minimized against {forcefield}")
        conformer_energies = mol.properties["conformer energies (kcal/mol)"]
        for conformer_idx, conformer_energy in enumerate(conformer_energies):
            print(
                "%5d / %5d : %8.3f kcal/mol %8.3f kcal/mol  %8.3f Angstroms"
                % (
                    conformer_idx + 1,
                    mol.n_conformers,
                    0.0,  # original energy goes here
                    conformer_energy / unit.kilocalories_per_mole,
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

    get_conformer_energies(
        molecule=args.molecule,
        registry=registry,
        forcefield=args.forcefield,
        constrained=args.constrained,
        prefix=args.prefix,
    )
