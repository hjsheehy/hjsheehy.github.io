import sys
import os
'''Generalise importing to modules with non .py extensions, e.g. .conf files.'''
from importlib.util import spec_from_loader, module_from_spec
from importlib.machinery import SourceFileLoader

def import_path(path):
    '''e.g. path = 'subfolder/file.conf' 
    Imports entire script with arbitrary file extension'''
    module_name = os.path.basename(path)#.replace('-', '_')
    spec = spec_from_loader(
        module_name,
        SourceFileLoader(module_name, path)
    )
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[module_name] = module
    return module
#
def import_all(module):
    '''Equivalent to 'from module import *' for modules with extensions
    other than .py'''
    module_dict = module.__dict__
    try:
        to_import = module.__all__
    except AttributeError:
        to_import = [name for name in module_dict if not name.startswith('_')]
    return module_dict, to_import
