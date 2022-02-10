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
def CreateName(model,hubbard,n,n_z,n_spins,n_orbs):
    return f'{model}_{hubbard}_{str(n).zfill(2)}_{str(n_z).zfill(2)}_{n_spins:d}_{n_orbs:d}'.replace('0.00e+00','0')


def Main():
    confname = os.path.join(out_folder,CreateName(model,hubbard,n,n_z,n_spins,n_orbs))
    temp_conf = f'''from lib import *
####################################################################
############################# Square ###############################
####################################################################
model='{model}'
n_name='{n_name}'
dos='True'
mu={mu}
t={t}
orbitals=[[0,0,0]]
n_orbs=len(orbitals)
n_spins={n_spins}
n_x=n_y={n}
n_z={n_z}
dimensions=[n_x,n_y,n_z]
V={V}
Delta_initial=1.42
phi_initial=0.97
hubbard='{hubbard}'
U={U}
T=0
friction={friction}
max_iterations=100
eps=0.001
'''
    with open(confname+'.conf', 'w') as f:
        f.write(temp_conf)
###################################################################
############################# Main ################################
###################################################################
n = 15
mu= -2.57
t = 1.0
n_orbs=1
V = 1.21
impurity_locations = [[0,0,0]]
epsilon = 0.1
U=2.8
friction=0.7

model_list=['tb','sc','sc']
spins_list=[1,2,2]
hubbard_list=['sc','nn','onsite']
n_name_list=['n','n_z']
n_z_list=[1,3]

ClearOutFolder(out_folder, silent)
for i,n_name in enumerate(n_name_list):
    n_z=n_z_list[i]
    for i in range(3):
        model=model_list[i]
        n_spins=spins_list[i]
        hubbard=hubbard_list[i]
        Main()
