
########################################################
######################### Import #######################
########################################################
import sys
import os
FILENAME=os.path.splitext(os.path.basename(sys.argv[0]))[0]
FILENAME=FILENAME.split('-')[0:2]
FILENAME=str('-'.join(FILENAME))

ROOT_DIR = os.path.dirname(
             os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT_DIR)
sys.path.append('01-main/')
from lib import *

silent=sys.argv[-1]

conf = os.path.join(ROOT_DIR,CONF,FILENAME+'*',OUT)
conf = glob.glob(conf)[0]

def conf_file(FILENAME):
    FILENAME_ = os.path.basename(FILENAME+'*')
    FILENAME_ = FILENAME_.split('.npz')[0]
    confname = os.path.join(conf,FILENAME_)+'.conf'
    config_module = import_path(os.path.join(MAIN,confname))
    module_dict, to_import = import_all(config_module)
    return {name: module_dict[name] for name in to_import}

########################################################
######################### Saving #######################
########################################################
OUTPUT = os.path.join(ROOT_DIR,FIG,FILENAME)

FILENAMES = glob.glob(DATA+FILENAME+'*/*')

def show_then_save(main, caption):
    fig=main()
    plt.show()

    fig=main()
    if silent!='-s':
        if YesNo('Save figure?'):
            plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight")
            caption()

########################################################
################### Phase diagram (2D) #################
########################################################
def name(filename):
    """Returns the filename without extension"""
    return os.path.basename(filename).split('.npz')[0]
def filename_x(filename):
    """Returns the x-coordinate as a float"""
    return float(filename.split('_')[2])
def filename_y(filename):
    """Returns the y-coordinate as a float"""
    return float(filename.split('_')[3])
def filename_z(filename):
    """Returns the z-coordinate, which is by default, the free energy"""
    return float(name(filename).split('_')[4])
def coincident_files(filename):
    """Returns all files which share the x and y coordinates of the filename"""
    filename=filename.split('/')[-1]
    x=filename_x(filename)
    y=filename_y(filename)
    filename=f'*_{x:.3f}_{y:.3f}_*'
    filenames = glob.glob(os.path.join(DATA,FILENAME,filename))
    return filenames
def minimum_z(filenames):
    """Returns the z-coordinate of the filename with the minimum z-coordinate"""
    free_energies=[filename_z(filename) for filename in filenames]
    return min(free_energies)
def filename_min_z(filenames):
    """Returns the filename of the filename with the minimum z-coordinate"""
    f=str(minimum_z(filenames))
    filename=f'*_{f}_*'
    filename = glob.glob(os.path.join(DATA,FILENAME,filename))
    filename=filename[0]
    return filename
def filenames_min_z():
    """Takes all FILENAMES and chooses the set with unique x, y and minimum z"""
    filenames=[]
    i=0
    for filename in FILENAMES:
        i+=1
        _filenames=coincident_files(filename)
        _filename=filename_min_z(_filenames)
        filenames=filenames+[_filename]
    filenames=list(set(filenames))
    return filenames
def sort_filenames(filenames):
    """Sorts the filenames row, column, in order to create a matrix"""
    return sorted(filenames, key=lambda k: [filename_y(k), filename_x(k)])
def filename_axes(filenames):
    """Finds the linspace lists of x and y coordinates"""
    x=sorted(set([filename_x(filename) for filename in filenames]))
    y=sorted(set([filename_y(filename) for filename in filenames]))
    return x,y
def matrix_filenames(filenames=FILENAMES):
    """Creates a matrix from the filenames"""
    xx,yy=filename_axes(filenames)
    lx=len(xx)
    ly=len(yy)
    col=[]
    for x in xx:
        row=[]
        for y in yy:
            filename=f'*_{x:.3f}_{y:.3f}_*'
            coincident = glob.glob(os.path.join(DATA,FILENAME+'*',filename))
            z=minimum_z(coincident)
            filename=f'*{x:.3f}_{y:.3f}_{z}*'

            filename = glob.glob(os.path.join(DATA,FILENAME+'*',filename))[0]
            row.append(filename)
        col.append(row)
    filenames=col
    return filenames

########################################################
########### Perculation of the phase diagram ###########
########################################################
def filename_seed(filename):
    """Returns the seed as an integer"""
    return int(name(filename).split('_')[1])
def filename_layer(filename):
    """Returns the layer as an integer"""
    return int(name(filename).split('_')[0])
