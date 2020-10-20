def _enforce_dependency_version(package, minimum_version):
    """
    Raise an exception if the toolkit version is less than a minimum version
    """
    import pkg_resources

    package_version = pkg_resources.get_distribution(package).version
    major, minor, patch = package_version.split(".")[:3]
    patch = patch.split("+")[0]
    minimum_version = minimum_version.split(".")[:3]
    assert all(
        [
            int(major) >= int(minimum_version[0]),
            int(minor) >= int(minimum_version[1]),
            int(patch) >= int(minimum_version[2]),
        ]
    ), f"Need at least version {'.'.join(minimum_version)} of package {package}"
