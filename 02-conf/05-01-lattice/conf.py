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
def CreateName(n):
        return f'lattice_{n:d}'

def Main():
    ClearOutFolder(out_folder)
    i=0
    confname = os.path.join(out_folder,CreateName(n))
    temp_conf = f'''n={n}
mu=-3.32
V=-1.21
v=w=1
zeta_v=0.52
zeta_w=0.23
Delta_v=0.92
Delta_w=0.73'''
    with open(confname+'.conf', 'w') as f:
        f.write(temp_conf)
###################################################################
############################# Main ################################
###################################################################
n=3
Main()
