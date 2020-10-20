import pytest


def test_dependency_version_check():
    """Ensure that runtime dependency checks are enforced"""
    from openff.cli.utils.utils import _enforce_dependency_version

    _enforce_dependency_version("openforcefield", "0.6.0")

    with pytest.raises(
        AssertionError, match="Need at least version 400.0.0 of package openforcefield"
    ):
        _enforce_dependency_version("openforcefield", "400.0.0")
