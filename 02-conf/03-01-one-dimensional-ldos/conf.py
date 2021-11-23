import sys
import os
ROOT_DIR = os.path.dirname(
            os.path.dirname(
             os.path.dirname(os.path.abspath(__file__))))
os.chdir(ROOT_DIR)
sys.path.append('01-main/')
from lib import *
fileName = os.path.dirname(os.path.realpath(__file__)).split('/')[-1]
out_folder = os.path.join(ROOT_DIR,CONF,fileName,OUT)
####################################################################
########################## Configuration ############################
####################################################################
def CreateName(n_x,phase):
    return f'Normal_{n_x:d}_{phase}'.replace('0.00e+00','0')

def Main():
    ClearOutFolder(out_folder)
    i=0
    for phase in phases:
        confname = os.path.join(out_folder,CreateName(n_x,phase))
        temp_conf = f'''from lib import *
####################################################################
######################## Self-consistent ###########################
####################################################################
model='tb'
lattice='SSH_intraorbital_nn'
phase='{phase}'
mu=0
n_spins=1
n_x={n_x}
V=0'''
        with open(confname+'.conf', 'w') as f:
            f.write(temp_conf)
        i+=1
    print(f'{i} .conf files created')
###################################################################
############################# Main ################################
###################################################################
n_x=15
phases=['lattice','trivial','trivial_t','trivial_tt','trivial_ttt','topological','topological_t','topological_tt','topological_ttt']
Main()
