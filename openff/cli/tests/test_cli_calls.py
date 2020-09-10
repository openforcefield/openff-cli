import subprocess


class TestCLICalls:
    def call(self, cmd):
        """Helper function to execute direct CLI calls"""
        if type(cmd) == str:
            cmd = cmd.split()

        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()

        return out.decode(), err.decode()

    def test_check_versions(self):
        """Test some basic behavior of check_versions.py"""
        from openforcefield import __version__ as toolkit_version

        from openff.cli import __version__ as cli_version

        out, _ = self.call("python openff/cli/check_versions.py")
        # TODO: Use regex to connect package names with versions
        assert toolkit_version in out
        assert cli_version in out
