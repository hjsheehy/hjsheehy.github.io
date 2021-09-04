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
def CreateName(model,n_name,n,n_z):
    return f'{model}_{n_name}_{str(n).zfill(2)}_{str(n_z).zfill(2)}'.replace('0.00e+00','0')


def Main():
    confname = os.path.join(out_folder,CreateName(model,n_name,n,n_z))
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
Delta_initial=0.45
phi_initial=2.085
hubbard='{hubbard}'
U={U}
T=0
friction=0
max_iterations=100
eps=0.001
'''
    with open(confname+'.conf', 'w') as f:
        f.write(temp_conf)
###################################################################
############################# Main ################################
###################################################################
n_list = np.arange(3,45,2)
n_z_list = np.arange(1,5,1)
mu= 3.57
t = 1.0
n_orbs=1
V = 0
impurity_locations = [[0,0,0]]
epsilon = 0.1
U=2.8

model_list=['tb','sc']
spins_list=[1,2]
hubbard='onsite'

ClearOutFolder(out_folder, silent)
for i in range(2):
    model=model_list[i]
    n_spins=spins_list[i]
    n_name='n'
    n_z=1
    for n in n_list:
        Main()
    n_name='n_z'
    n=41
    for n_z in n_z_list:
        Main()
