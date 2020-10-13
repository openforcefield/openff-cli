import pathlib
import subprocess
import tempfile

from openforcefield.utils import temporary_cd

from openff.cli import __file__ as cli_path

# TODO: Run tests from tempdirs
CLI_ROOT = pathlib.Path(cli_path).joinpath("../../../").resolve()


class TestCLICalls:
    def call(self, cmd, raise_err=True):
        """Helper function to execute direct CLI calls"""
        if type(cmd) == str:
            cmd = cmd.split()

        with tempfile.TemporaryDirectory() as tmp_dir:
            with temporary_cd(tmp_dir):
                proc = subprocess.Popen(
                    cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                out, err = proc.communicate()

        if err and raise_err:
            raise Exception(err)
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
            f"python {CLI_ROOT}/openff/cli/generate_conformers.py --toolkit magic",
            raise_err=False,
        )
        assert "openff.cli.utils.exceptions.UnsupportedToolkitError" in err
        assert "magic" in err
