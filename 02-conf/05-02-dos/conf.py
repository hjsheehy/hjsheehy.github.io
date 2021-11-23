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
mu=3.92
V=2.21
w=0.79
v=1.31
w=v=1
U_T=2.56
U_R=0.12
phi=2.32
zeta_v=2.52
zeta_w=1.23
Delta_v=3.92
Delta_w=2.73'''
    with open(confname+'.conf', 'w') as f:
        f.write(temp_conf)
###################################################################
############################# Main ################################
###################################################################
n=41
Main()
