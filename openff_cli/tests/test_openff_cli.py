"""
Unit and regression test for the openff_cli package.
"""

# Import package, test suite, and other packages as needed
import openff_cli
import pytest
import sys

def test_openff_cli_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "openff_cli" in sys.modules
