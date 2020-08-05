import sys

from openff.cli import check_versions


class TestCheckVersions:
    def test_basic_output(self):
        out = check_versions.get_versions()

        assert "OpenFF Toolkit" in out
        assert "OpenFF CLI" in out
