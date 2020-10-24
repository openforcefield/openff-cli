def _enforce_dependency_version(package, minimum_version):
    """
    Raise an exception if the toolkit version is less than a minimum version
    """
    import pkg_resources
    from packaging.version import parse

    package_version = pkg_resources.get_distribution(package).version

    assert parse(package_version) > parse(
        minimum_version
    ), "Need at least version {'.'.join(minimum_version)} of package {package}"
