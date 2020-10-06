import pathlib
import tempfile

import numpy as np
import pytest
from openforcefield.tests.utils import requires_openeye, requires_rdkit
from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff.forcefield import get_available_force_fields
from openforcefield.utils import (
    OpenEyeToolkitWrapper,
    RDKitToolkitWrapper,
    temporary_cd,
)

from openff.cli.check_versions import get_versions
from openff.cli.generate_conformers import (
    generate_conformers,
    make_registry,
    write_mols,
)
from openff.cli.tests.utils import get_data_file_path

# TODO: Run all tests in a safe temporary directory


class TestCheckVersions:
    def test_basic_output(self):
        out = get_versions()
        from openforcefield import __version__ as toolkit_version

        from openff.cli import __version__ as cli_version

        assert "OpenFF Toolkit" in out
        assert "OpenFF CLI" in out
        assert toolkit_version in out
        assert cli_version in out


@pytest.mark.parametrize(
    "toolkit",
    [
        pytest.param("rdkit", marks=requires_rdkit),
        pytest.param("openeye", marks=requires_openeye),
    ],
)
class TestGenerateConformersCLI:
    def test_load_one_mol_sdf_without_charge(self, toolkit):
        """Test loading one molecule from a .sdf file WITHOUT charges"""
        registry = make_registry(toolkit)
        ethanol = get_data_file_path("molecules/ethanol.sdf")
        assert Molecule.from_file(ethanol).partial_charges is None

        mols_out = generate_conformers(
            molecule=ethanol,
            forcefield="openff-1.0.0.offxml",
            registry=registry,
        )

        assert len(mols_out) > 1
        assert mols_out[0].partial_charges is not None

    def test_load_one_mol_sdf_with_charge(self, toolkit):
        """Test loading one molecule from a .sdf file WITH charges"""
        registry = make_registry(toolkit)
        ethanol_partial_charges = get_data_file_path(
            "molecules/ethanol_partial_charges.sdf"
        )
        charges_in = Molecule.from_file(ethanol_partial_charges).partial_charges
        assert charges_in is not None

        mols_out = generate_conformers(
            molecule=ethanol_partial_charges,
            forcefield="openff-1.0.0.offxml",
            registry=registry,
        )

        assert len(mols_out) > 1
        assert np.allclose(mols_out[0].partial_charges, charges_in)

    # TODO: Figure out how to skip a test based on a variable passed in through parametrize
    @pytest.mark.skipif(True, reason="Test requires OpenEye toolkit")
    def test_load_one_mol_mol2_without_charge(self, toolkit):
        """Test loading one molecule from a .mol2 file WITH charges"""
        registry = make_registry(toolkit)
        toluene_partial_charges = get_data_file_path("molecules/toluene_charged.mol2")
        charges_in = Molecule.from_file(toluene_partial_charges).partial_charges

        mols_out = generate_conformers(
            molecule=toluene_partial_charges,
            forcefield="openff-1.0.0.offxml",
            registry=registry,
        )

        assert len(mols_out) > 1
        assert np.allclose(mols_out[0].partial_charges, charges_in)

    @pytest.mark.skipif(True, reason="Test requires OpenEye toolkit")
    def test_load_one_mol_mol2_with_charge(self, toolkit):
        """
        Test loading one molecule from a .mol2 file WITHOUT charges

        .. note :: This file was generated via the one-liner below and has nan
            charges, which may or may not be valid.

        ```
        Molecule.from_smiles('CCO').to_file(
            'ethanol_no_charges.mol2',
            file_format='mol2',
        )
        ```
        """
        registry = make_registry(toolkit)
        ethanol_partial_charges = get_data_file_path(
            "molecules/ethanol_no_charges.mol2"
        )
        charges_in = Molecule.from_file(ethanol_partial_charges).partial_charges

        mols_out = generate_conformers(
            molecule=ethanol_partial_charges,
            forcefield="openff-1.0.0.offxml",
            registry=registry,
        )

        assert len(mols_out) == 1

        assert not charges_in
        assert not mols_out[0].partial_charges

    # loading a molecule with a defined name (from any format) and ensuring the output file has that prefix
    # loading a molecule with a defined name (from any format), and providing
    # a -f option and ensuring the output file has the -f prefix

    def test_p_prefix(self, toolkit):
        """Ensure the output file has the -p prefix"""
        registry = make_registry(toolkit)

        ethanol = get_data_file_path("molecules/ethanol.sdf")

        mols_out = generate_conformers(
            molecule=ethanol,
            forcefield="openff-1.0.0.offxml",
            registry=registry,
            prefix="test_ethanol",
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            with temporary_cd(tmp_dir):
                write_mols(
                    mols=mols_out,
                    toolkit_registry=registry,
                    molecule=ethanol,
                    prefix="test_ethanol",
                )

                assert pathlib.Path("test_ethanol_0.sdf").is_file()

    # loading a molecule with a defined name (from any format) and ensuring the output file has that prefix

    def test_load_multi_mol_sdf(self, toolkit):
        """Test the case of an SDF file with multiple molecules"""
        registry = make_registry(toolkit)
        butane_multi = get_data_file_path("molecules/butane_multi.sdf")
        generate_conformers(
            molecule=butane_multi,
            forcefield="openff-1.0.0.offxml",
            registry=registry,
        )

    def test_load_one_mol_smi(self, toolkit):
        """Test loading one molecule from SMILES in a .smi file"""
        registry = make_registry(toolkit)
        ebastine = get_data_file_path("molecules/ebastine.smi")
        mols_out = generate_conformers(
            molecule=ebastine,
            forcefield="openff-1.0.0.offxml",
            registry=registry,
        )

        assert len(mols_out) > 1

    def test_load_multi_mol_smi(self, toolkit):
        """Test loading multiple molecules from SMILES in a .smi file"""
        registry = make_registry(toolkit)
        dyes = get_data_file_path("molecules/multi_mols.smi")
        mols_out = generate_conformers(
            molecule=dyes,
            forcefield="openff-1.0.0.offxml",
            registry=registry,
        )

        assert len(mols_out) > 3

    def test_load_ambiguous_stereo_smi(self, toolkit):
        """Test loading a molecule with ambiguous stereo from SMILES and enumerating stereoisomers"""
        # TODO: Should the CLI accept both paths and SMILES as strings, or only files?
        registry = make_registry(toolkit)
        mol = get_data_file_path("molecules/dichloroethene_ambiguous_stereo.smi")
        mols_out = generate_conformers(
            molecule=mol,
            forcefield="openff-1.0.0.offxml",
            registry=registry,
        )

        assert len(mols_out) == 2

    def test_preserve_stereo_smi(self, toolkit):
        """Test loading a molecule with defined stereo from SMILES and preserving that stereochemistry"""
        registry = make_registry(toolkit)
        mol = get_data_file_path("molecules/dichloroethene_stereo.smi")
        mols_out = generate_conformers(
            molecule=mol,
            forcefield="openff-1.0.0.offxml",
            registry=registry,
        )

        assert Molecule(mol).is_isomorphic_with(mols_out[0])

    @pytest.mark.parametrize("ff_name", get_available_force_fields())
    def test_available_forcefields(self, toolkit, ff_name):
        """Test that all forcefields available through get_available_forcefields
        can at least do a basic conformer generaiton"""
        # TODO: Possibly ensure this produces non-garbage output,
        #  instead of only "doesn't raise an error we know of"
        registry = make_registry(toolkit)
        mol = get_data_file_path("molecules/ethanol.sdf")
        generate_conformers(
            molecule=mol,
            forcefield=ff_name,
            registry=registry,
        )

    def test_different_forcefields(self, toolkit):
        """Test that different --forcefield arguments produce different results"""
        from openff.cli.tests.utils import get_data_file_path

        registry = make_registry(toolkit)
        mol = get_data_file_path("molecules/ebastine.smi")
        parsley_1_0_0 = generate_conformers(
            molecule=mol,
            forcefield="openff-1.0.0",
            prefix="parsley100",
            registry=registry,
        )
        parsley_1_2_0 = generate_conformers(
            molecule=mol,
            forcefield="openff-1.2.0",
            prefix="parsley120",
            registry=registry,
        )

        assert not np.allclose(
            parsley_1_0_0[0].conformers[0],
            parsley_1_2_0[0].conformers[0],
        )


@requires_rdkit
@requires_openeye
def test_make_registry(self, toolkit):
    """Test the behavior of the make_registry helper function. This is
    hard-coded to each toolkit as specific objects needed to be found."""
    rdkit_registry = make_registry("rdkit")
    assert isinstance(rdkit_registry.registered_toolkits[0], RDKitToolkitWrapper)
    assert OpenEyeToolkitWrapper not in rdkit_registry.registered_toolkits

    openeye_registry = make_registry("openeye")
    assert isinstance(openeye_registry.registered_toolkits[0], OpenEyeToolkitWrapper)
    assert RDKitToolkitWrapper not in rdkit_registry.registered_toolkits


@requires_rdkit
@requires_openeye
def test_different_toolkits():
    """Test that different --toolkit argumnents produce different results"""
    from openff.cli.tests.utils import get_data_file_path

    rdkit = make_registry("rdkit")
    openeye = make_registry("openeye")
    mol = get_data_file_path("molecules/ebastine.smi")
    rdkit_results = generate_conformers(
        molecule=mol,
        forcefield="openff-1.0.0",
        prefix="rdkit",
        registry=rdkit,
    )
    openeye_results = generate_conformers(
        molecule=mol,
        forcefield="openff-1.0.0",
        prefix="openeye",
        registry=openeye,
    )

    assert not np.allclose(
        rdkit_results[0].conformers[0],
        openeye_results[0].conformers[0],
    )


def test_parsing_error():
    """Ensure that failed molecule parsing is caught"""
    from openff.cli.utils.exceptions import MoleculeParsingError

    # Ensure the exception can be raised with no keyword arguments
    exec_msg = "Failed to parse a molecule file."
    with pytest.raises(MoleculeParsingError, match=exec_msg):
        raise MoleculeParsingError

    # Try to parse a MOL2 file with RDKit
    registry = make_registry("rdkit")
    mol = get_data_file_path("molecules/ethanol_no_charges.mol2")
    exec_msg = (
        r"Failed to parse a molecule file. Attempted to parse file "
        ".*ethanol_no_charges.mol2 using toolkit registry containing+.*RDKit"
    )
    with pytest.raises(MoleculeParsingError, match=exec_msg):
        generate_conformers(
            molecule=mol,
            forcefield="openff-1.0.0.offxml",
            registry=registry,
        )


def test_unsupported_toolkit_error():
    """Ensure that requesting an unsupported toolkit is caught"""
    from openff.cli.utils.exceptions import UnsupportedToolkitError

    exec_msg = (
        r"Requested use of a cheminformatics toolkit not supported by "
        "the OpenFF Toolkit. .* magic"
    )
    with pytest.raises(UnsupportedToolkitError, match=exec_msg):
        make_registry("magic")
