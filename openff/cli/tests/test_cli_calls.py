import pathlib
import subprocess
import tempfile

from openforcefield.utils import temporary_cd

from openff.cli import __file__ as cli_path

CLI_ROOT = pathlib.Path(cli_path).joinpath("../../../").resolve()


class TestCLICalls:
    def call(self, cmd):
        """Helper function to execute direct CLI calls"""
        if type(cmd) == str:
            cmd = cmd.split()

        with tempfile.TemporaryDirectory() as tmp_dir:
            with temporary_cd(tmp_dir):
                proc = subprocess.Popen(
                    cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                out, err = proc.communicate()

        return out.decode(), err.decode()


class TestCheckVersionsCalls(TestCLICalls):
    def test_check_versions(self):
        """Test some basic behavior of check_versions.py"""
        from openforcefield import __version__ as toolkit_version

        from openff.cli import __version__ as cli_version

        out, _ = self.call(f"python {CLI_ROOT}/openff/cli/check_versions.py")
        # TODO: Use regex to connect package names with versions
        assert toolkit_version in out
        assert cli_version in out


class TestGenerateConformersCalls(TestCLICalls):
    def test_unsupported_toolkit(self):
        """Ensure that requesting an unsupported toolkit is caught"""
        # TODO: This should maybe fail before checking the toolkit,
        #  based on missing required arguments
        _, err = self.call(
            f"python {CLI_ROOT}/openff/cli/generate_conformers.py --toolkit magic"
        )
        assert "openff.cli.utils.exceptions.UnsupportedToolkitError" in err
        assert "magic" in err

    def test_different_forcefields(self):
        """Test that different --forcefield arguments produce different results"""
        import numpy as np
        from openforcefield.topology.molecule import Molecule

        from openff.cli.tests.utils import get_data_file_path

        filepath = get_data_file_path("molecules/ethanol.sdf")
        cmd = f"python {CLI_ROOT}/openff/cli/generate_conformers.py --toolkit rdkit"
        cmd += f" --constrained False --molecule {filepath}"
        parsley_1_0_0 = cmd + " --forcefield openff-1.0.0 --prefix parsley120"
        parsley_1_2_0 = cmd + " --forcefield openff-1.2.0 --prefix parsley100"

        self.call(parsley_1_0_0)
        self.call(parsley_1_2_0)

        assert not np.allclose(
            Molecule.from_file("parsley120_0.sdf").conformers[0],
            Molecule.from_file("parsley100_0.sdf").conformers[0],
        )

    def test_different_toolkits(self):
        """Test that different --toolkit arguments produce different results"""
        import numpy as np
        from openforcefield.topology.molecule import Molecule

        from openff.cli.tests.utils import get_data_file_path

        # Use a large molecule to ensure conformer generation results differ
        filepath = get_data_file_path("molecules/ebastine.smi")
        cmd = f"python {CLI_ROOT}/openff/cli/generate_conformers.py --forcefield openff-1.0.0"
        cmd += f" --constrained False --molecule {filepath}"
        rdkit = cmd + " --toolkit rdkit --prefix rdkit"
        openeye = cmd + " --toolkit openeye --prefix openeye"

        self.call(rdkit)
        self.call(openeye)

        assert not np.allclose(
            Molecule.from_file("rdkit_0.sdf").conformers[0],
            Molecule.from_file("openeye_0.sdf").conformers[0],
        )
