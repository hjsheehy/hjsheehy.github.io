#########################################################
################# from .conf import * ###################
#########################################################
from lib import *
from phase_diagram import *

config_file=True

"""Imports .conf file. Must be ran as "python example_simulation.py example_conf.py" """
global CONFNAME, SIM_NAME

if len(sys.argv) < 2:
        print('Please supply conf file')
        sys.exit()
    
conf = sys.argv[-1]

CONFNAME = conf.split('.conf')[0]
config_module = import_path(os.path.join(MAIN,conf))
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

CONFNAME = os.path.basename(CONFNAME) 
SIM_NAME = os.path.dirname(sys.argv[1])
SIM_NAME = os.path.basename(SIM_NAME)
CONFNAME = os.path.join(DATA,SIM_NAME,CONFNAME)
