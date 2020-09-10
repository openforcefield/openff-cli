class MoleculeParsingError(Exception):
    """An exception raised when the OpenFF Toolkit appears to fail while parsing molecule files."""

    def __init__(self, filename=None, toolkit_registry=None):
        """
        Parameters
        ----------
        toolkit_registry: ToolkitRegistry
            The registry which attempted to parse the molecule
            source.
        """
        self.message = "Failed to parse a molecule file."
        if filename is not None:
            self.message += f" Attempted to parse file {filename}"
        if toolkit_registry is not None:
            self.message += (
                " using toolkit registry containing the following toolkits: "
            )
            # This attrib was only added in toolkit > 0.7.1
            if hasattr(toolkit_registry, "registered_toolkit_versions"):
                self.message += str(toolkit_registry.registered_toolkit_versions)
            else:
                self.message += str(toolkit_registry.registered_toolkits)

    def __str__(self):
        return self.message
