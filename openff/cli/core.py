from typing import Tuple

import numpy as np
from openforcefield.topology.molecule import Molecule
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.utils.toolkits import ToolkitRegistry
from simtk import openmm, unit


def make_registry(toolkit: str) -> ToolkitRegistry:
    if toolkit.lower() == "openeye":
        from openforcefield.utils.toolkits import OpenEyeToolkitWrapper

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
    elif toolkit.lower() == "rdkit":
        from openforcefield.utils.toolkits import RDKitToolkitWrapper

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper])
    else:
        from openff.cli.utils.exceptions import UnsupportedToolkitError

        raise UnsupportedToolkitError(toolkit=toolkit)

    # Checks later assume that this is length 1. This should be changed if
    # multiple toolkits (i.e. RDKit and AmberTools) are needed at once
    assert len(toolkit_registry.registered_toolkit_versions) == 1
    return toolkit_registry


def _get_forcefield(forcefield: str, constrained: bool) -> ForceField:
    ff_name = forcefield
    if constrained:
        split_name = ff_name.split("-")
        split_name[0] = "openff_unconstrained"
        ff_name = "-".join(split_name)
    if not ff_name.endswith(".offxml"):
        ff_name += ".offxml"
    ff = ForceField(ff_name)

    return ff


def _build_simulation(molecule, forcefield, mols_with_charge=None):
    """Given a Molecule and ForceField, initialize a barebones OpenMM Simulation."""
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
            for i in range(molecule.n_atoms)
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


def _get_conformer_data(
    simulation: openmm.app.Simulation,
) -> Tuple[unit.Quantity, unit.Quantity]:
    """Given an OpenMM simulation, return its energy and coordinates"""
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy()
    coords = state.getPositions()

    coords = np.array([[atom.x, atom.y, atom.z] for atom in coords]) * unit.nanometer

    return energy, coords


def _minimize_conformer(
    simulation: openmm.app.Simulation, conformer: unit.Quantity
) -> openmm.app.Simulation:
    """Minimize a conformer via OpenMM"""
    simulation.context.setPositions(conformer)
    simulation.minimizeEnergy()

    return simulation


def _get_rms_two_conformers(
    mol: Molecule, positions1: unit.Quantity, positions2: unit.Quantity
) -> float:
    """Find the RMSD between two conformers of a molecule using RDKit"""
    # TODO: Is it worth making Molecule.get_rmsd(), which operates
    # through ToolkitWrapper methods?
    from rdkit.Chem import rdMolAlign

    mol_copy = Molecule(mol)
    mol_copy._conformers = None
    mol_copy.add_conformer(positions1)
    mol_copy.add_conformer(positions2)

    rdmol = mol_copy.to_rdkit()
    rmslist = []
    rdMolAlign.AlignMolConformers(rdmol, RMSlist=rmslist)

    return rmslist[0]
