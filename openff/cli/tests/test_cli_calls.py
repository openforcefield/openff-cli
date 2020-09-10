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
        import openforcefield

        import openff.cli

        out, _ = self.call("python openff/cli/check_versions.py")
        # TODO: Use regex to connect package names with versions
        assert openforcefield.__version__ in out
        assert openff.cli.__version__ in out
