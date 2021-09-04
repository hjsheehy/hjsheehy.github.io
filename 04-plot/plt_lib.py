########################################################
######################### Import #######################
########################################################
import sys
import os
FILENAME=os.path.splitext(os.path.basename(sys.argv[0]))[0]

ROOT_DIR = os.path.dirname(
             os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT_DIR)
sys.path.append('01-main/')
from lib import *

silent=sys.argv[-1]

conf = os.path.join(ROOT_DIR,CONF,FILENAME,OUT)

def conf_file(filename):
    filename_ = os.path.basename(filename).split('.npz')[0]
    confname = os.path.join(conf,filename_)+'.conf'
    config_module = import_path(os.path.join(MAIN,confname))
    module_dict, to_import = import_all(config_module)
    return {name: module_dict[name] for name in to_import}

########################################################
######################### Saving #######################
########################################################
output = os.path.join(ROOT_DIR,FIG,FILENAME)

filenames = glob.glob(DATA+FILENAME+'/*')

def show_then_save(main, caption):
    fig=main()
    plt.show()

    fig=main()
    if silent!='-s':
        if YesNo('Save figure?'):
            plt.savefig(output+'.pdf', bbox_inches = "tight")
            caption()

########################################################
################## Plotting functions ##################
########################################################

########################################################
######################### Config #######################
########################################################
latex_width=4.7747
