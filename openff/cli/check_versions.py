import importlib


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

    print(out)


if __name__ == "__main__":
    # TODO: Have this function return the string, without printing above
    #  or return an exit code here
    get_versions()
