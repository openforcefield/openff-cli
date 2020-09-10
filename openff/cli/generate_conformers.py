import argparse
from copy import deepcopy
from typing import List, Optional, Tuple

# from rdkit.Chem import rdMolAlign
import numpy as np
from openforcefield.topology.molecule import Molecule, UndefinedStereochemistryError
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.utils.toolkits import ToolkitRegistry
from simtk import openmm, unit


def generate_conformers(
    molecule: str,
    registry: ToolkitRegistry,
    forcefield: str,
    constrained: bool = False,
    prefix: Optional[str] = None,
) -> List[Molecule]:

    ff_name = forcefield
    if constrained:
        split_name = ff_name.split("-")
        split_name[0] = "openff_unconstrained"
        ff_name = "-".join(split_name)
    if not ff_name.endswith(".offxml"):
        ff_name += ".offxml"
    ff = ForceField(ff_name)

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

    # If no conformers were found (i.e. SMI files), generate some
    # TODO: How many should be generated?
    # TODO: If 1 or a few conformers are found, should more be generated?
    for mol in mols:
        if mol.conformers is None:
            mol.generate_conformers(toolkit_registry=registry, n_conformers=1)

    # TODO: What happens if some molecules in a multi-molecule file have charges, others don't?
    mols_with_charges = []
    for mol in mols:
        if mol.partial_charges is not None:
            mols_with_charges.append(mol)

    mols = _collapse_conformers(mols)

    mols_out = []
    for mol in mols:
        simulation, partial_charges = _build_simulation(
            molecule=mol, forcefield=ff, mols_with_charge=mols_with_charges
        )
        mol._partial_charges = partial_charges

        for i, conformer in enumerate(mol.conformers):
            energy, positions = _get_minimized_data(conformer, simulation)
            mol = _reconstruct_mol_from_conformer(mol, positions)
            _add_metadata_to_mol(mol, energy, registry.toolkit_version, ff_name)
            mol.name += "_conf" + str(i)
            mols_out.append(mol)

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


def _build_simulation(molecule, forcefield, mols_with_charge):
    """Given a Molecule and ForceField, initialize a barebones OpenMM Simulation."""
    mol_copy = deepcopy(molecule)
    off_top = molecule.to_topology()
    system, ret_top = forcefield.create_openmm_system(
        off_top,
        charge_from_molecules=mols_with_charge,
        allow_nonintegral_charges=True,
        return_topology=True,
    )

    # Use OpenMM to compute initial and minimized energy for all conformers
    integrator = openmm.VerletIntegrator(1 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName("Reference")
    omm_top = off_top.to_openmm()
    simulation = openmm.app.Simulation(omm_top, system, integrator, platform)

    charges_from_top = [*ret_top.reference_molecules][0].partial_charges
    if charges_from_top is not None:
        partial_charges = charges_from_top
    else:
        # ret_top only has partial charges in OFFTK 0.8.0+, so may need to
        # manually get partial charges from OpenMM if using OFFTK <= 0.7.1
        partial_charges = [
            system.getForces()[0].getParticleParameters(i)[0]
            for i in range(mol_copy.n_atoms)
        ]
        # Unwrap list of Quantity objects into a single Quantity that contains a list
        # Surely there's a simpler way to to this?
        partial_charges = unit.Quantity(
            np.asarray(
                [val.value_in_unit(unit.elementary_charge) for val in partial_charges]
            ),
            unit=unit.elementary_charge,
        )

    return simulation, partial_charges


def _get_minimized_data(
    conformer: unit.Quantity, simulation: openmm.app.Simulation
) -> Tuple[unit.Quantity]:
    """Given an OpenMM simulation and conformer, minimze and return an energy"""
    simulation.context.setPositions(conformer)
    simulation.minimizeEnergy()

    min_state = simulation.context.getState(getEnergy=True, getPositions=True)
    min_energy = min_state.getPotentialEnergy()
    min_coords = min_state.getPositions()

    return min_energy, min_coords


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
    mol: Molecule, energy: unit.Quantity, toolkit_version: str, ff_name: str
) -> None:
    mol.properties["absolute energy (kcal/mol): "] = energy
    mol.properties["conformer generation toolkit: "] = toolkit_version
    mol.properties["minimized against: "] = ff_name


def make_registry(toolkit: str) -> ToolkitRegistry:
    if toolkit.lower() == "openeye":
        import openeye
        from openforcefield.utils.toolkits import OpenEyeToolkitWrapper

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        toolkit_version = "openeye-toolkits version " + openeye.__version__
    elif toolkit.lower() == "rdkit":
        import rdkit
        from openforcefield.utils.toolkits import RDKitToolkitWrapper

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper])
        toolkit_version = "rdkit version " + rdkit.__version__
    else:
        from openff.cli.utils.exceptions import UnsupportedToolkitError

        raise UnsupportedToolkitError(toolkit=toolkit)

    if not hasattr(toolkit_registry, "toolkit_version"):
        toolkit_registry.toolkit_version = toolkit_version

    return toolkit_registry


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
    parser.add_argument(
        "-t",
        "--toolkit",
        type=str,
        help="Name of the underlying cheminformatics toolkit to use",
    )
    parser.add_argument(
        "-f",
        "--forcefield",
        type=str,
        help="Name of the force field to use",
    )
    parser.add_argument(
        "-m",
        "--molecule",
        type=str,
        help="Path to an input file containing a molecule(s)",
    )
    parser.add_argument(
        "-r",
        "--rms-cutoff",
        type=float,
        default=0.25,
        help="The redundancy cutoff between pre-minimized conformers",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        help="The prefix for output molecules",
    )
    parser.add_argument(
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
