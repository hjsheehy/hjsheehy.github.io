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
def CreateName(model,hubbard,n,n_spins,n_orbs):
    return f'{model}_{hubbard}_{str(n).zfill(2)}_{n_spins:d}_{n_orbs:d}'.replace('0.00e+00','0')


def Main():
    ClearOutFolder(out_folder, silent)
    for i in range(3):
        model=model_list[i]
        n_spins=spins_list[i]
        hubbard=hubbard_list[i]
        confname = os.path.join(out_folder,CreateName(model,hubbard,n,n_spins,n_orbs))
        temp_conf = f'''from lib import *
####################################################################
############################# Square ###############################
####################################################################
model='{model}'
dos='True'
mu={mu}
t={t}
orbitals=[[0,0]]
n_orbs=len(orbitals)
n_spins={n_spins}
n_x=n_y={n}
dimensions=[n_x,n_y]
V={V}
Delta_initial=1.42
phi_initial=0.97
hubbard='{hubbard}'
U={U}
T=0
friction=0.
max_iterations=100
eps=0.001
'''
        with open(confname+'.conf', 'w') as f:
            f.write(temp_conf)
###################################################################
############################# Main ################################
###################################################################
n = 53
mu= -3.57
t = 1.0
n_orbs=1
V = 1.21
impurity_locations = [[0,0,0]]
epsilon = 0.1
U=2.8

model_list=['tb','sc','sc']
spins_list=[1,2,2]
hubbard_list=['sc','nn','onsite']


Main()
