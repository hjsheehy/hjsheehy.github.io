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

silent=False
if sys.argv[-1]=='-s':
    silent=True
####################################################################
########################## Configuration ############################
####################################################################
def CreateName(hubbard,n,n_spins,n_orbs):
    return f'{hubbard}_{str(n).zfill(2)}_{n_spins:d}_{n_orbs:d}'.replace('0.00e+00','0')


def Main():
    ClearOutFolder(out_folder, silent)
    for i in range(len(hubbard_list)):
        model=model_list[i]
        hubbard=hubbard_list[i]
        confname = os.path.join(out_folder,CreateName(hubbard,n,n_spins,n_orbs))
        temp_conf = f'''from lib import *
####################################################################
############################# Square ###############################
####################################################################
model='{model}'
mu={mu}
t={t}
orbitals=[[0,0]]
n_orbs=len(orbitals)
n_spins={n_spins}
n_x=n_y={n}
n_z=1
V={V}
Delta=1.42
phi=.97
zeta=0.11
hubbard='{hubbard}'
U={U}
U_nn={U_nn}
U_nn_opp={U_nn_opp}
T=0
friction=0.4
max_iterations=100
eps=0.001
'''
        with open(confname+'.conf', 'w') as f:
            f.write(temp_conf)
###################################################################
############################# Main ################################
###################################################################
n = 53
mu= -3.27
t = 1.0
n_orbs=1
n_spins=2
V = 1.21
impurity_locations = [[0,0,0]]
epsilon = 0.1
U=2.8
U_nn=1.2
U_nn_opp=0.8

model_list=['tb','sc','sc','sc']
hubbard_list=['tb','nn_equal','nn_opposite','onsite']


Main()
