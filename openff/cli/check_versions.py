import importlib


def get_versions():
    packages = {
        'OpenFF Toolkit': 'openforcefield',
        'Openforcefields': 'openforcefields',
        'OpenFF Evaluator': 'openff.evaluator',
        'OpenFF System': 'openff.system',
        'OpenFF CLI': 'openff.cli',
        'CMILES': 'cmiles',
    }

    out = 'Found the following packages installed\n' + \
            '{0:20}\t{1}\n'.format('Package name', 'Version') + \
           '------------            -------\n'

    for key, val in packages.items():
        try:
            mod = importlib.import_module(val)
        except ImportError:
            out += f'{key:20}\t' 'Not found\n'
            continue
        out += f'{key:20}\t{mod.__version__}\n'

    return out

 
if __name__ == '__main__':
    exit(get_versions())
