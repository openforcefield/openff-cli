import importlib
import logging

import click

# The OpenFF toolkit's use of the logging modules prints its "you don't have
# OpenEye installed" error to STDERR. We do not wish to capture this here, so
# set the minimum logging level to ERROR. See
# https://docs.python.org/3/library/logging.html#logging-levels

logging.getLogger("openforcefield").setLevel(logging.ERROR)


def get_versions():
    packages = {
        "OpenFF Toolkit": "openforcefield",
        "Openforcefields": "openforcefields",
        "OpenFF Evaluator": "openff.evaluator",
        "OpenFF System": "openff.system",
        "OpenFF CLI": "openff.cli",
        "CMILES": "cmiles",
    }

    out = (
        "Found the following packages installed\n"
        "{0:20}\t{1}\n------------            -------\n".format(
            "Package name", "Version"
        )
    )

    for key, val in packages.items():
        try:
            mod = importlib.import_module(val)
        except ImportError:
            out += f"{key:20}\t" "Not found\n"
            continue
        out += f"{key:20}\t{mod.__version__}\n"

    return out


@click.command("check_versions")
def check_versions_command():
    """Check the installed versions of OpenFF software."""
    click.echo(get_versions())
