import pkg_resources


def get_data_file_path(rel_path):
    """
    Get the path to a data file. Normally these can be found in
    `openff/cli/data/` but are likely hidden somewhere else in site-packages
    after installation.

    Parameters
    ----------
    file_path : str
        Name of the file to find the full path of
    """

    full_path = pkg_resources.resource_filename("openff.cli", "data/" + rel_path)

    return full_path
