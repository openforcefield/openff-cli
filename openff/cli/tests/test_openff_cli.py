"""
Unit and regression test for the openff_cli package.
"""

# Import package, test suite, and other packages as needed
from openff import cli
import pytest
import sys

def test_openff_cli_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "openff.cli" in sys.modules
