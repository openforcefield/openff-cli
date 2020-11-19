import click

from openff.cli.check_versions import check_versions_command
from openff.cli.generate_conformers import generate_conformers_command
from openff.cli.get_conformer_energies import get_conformer_energies_command


@click.group()
def cli():
    """Command line utilities for the Open Force Field software stack"""


cli.add_command(check_versions_command)
cli.add_command(generate_conformers_command)
cli.add_command(get_conformer_energies_command)
