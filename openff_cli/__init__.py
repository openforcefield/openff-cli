"""
OpenFF CLI
Command line utilities for the Open Force Field software stack
"""

# Add imports here
from .cli import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
