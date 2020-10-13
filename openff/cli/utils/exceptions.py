class MoleculeParsingError(Exception):
    """An exception raised when the OpenFF Toolkit appears to fail while parsing molecule files."""

    def __init__(self, filename=None, toolkit_registry=None):
        """
        Parameters
        ----------
        filename : str, optional
            The name of the file that was requested to be loaded in
        toolkit_registry : ToolkitRegistry, optional
            The registry which attempted to parse the molecule source.
        """
        self.message = "Failed to parse a molecule file."
        if filename:
            self.message += f" Attempted to parse file {filename}"
        if toolkit_registry:
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


class UnsupportedToolkitError(Exception):
    """An exception raised when requesting use of a cheminformatics toolkit
    not supported by the OpenFF Toolkit."""

    def __init__(self, toolkit=None):
        """
        Parameters
        ----------
        toolkit : str, optional
            The name of the cheminformatics toolkit that was requested
        """
        self.message = (
            "Requested use of a cheminformatics toolkit not supported by "
            "the OpenFF Toolkit."
        )
        if toolkit:
            self.message += f" Requested the toolkit of name {toolkit}"

    def __str__(self):
        return self.message
