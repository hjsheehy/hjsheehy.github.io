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
def CreateName(n_x,n_y,n_z,mu,t,V):
        return f'Singlet_{n_x:d}_{n_y:d}_{n_z:d}_{mu:0.2f}_{t:0.2f}_{V:0.2e}'.replace('0.00e+00','0')

def Main():
    ClearOutFolder(out_folder, silent)
    for n in n_list:
        n_z=n    
        confname = os.path.join(out_folder,CreateName(n_x,n_y,n_z,mu,t,V))
        temp_conf = f'''from lib import *
####################################################################
############################# Square ###############################
####################################################################
model='bdg_singlet_sc'
n_spins={n_spins}
n_x={n_x}
n_y={n_y}
n_z={n_z}
a={a}

a1=a*np.array([1,0,0])
a2=a*np.array([0,1,0])
a3=a*np.array([0,0,1])
zero=np.array([0,0,0])
dimensions=[n_x,n_y,n_z]
basis=[a1,a2,a3]
######### lattice & sites #########
mu={mu}
orbitals=[zero]
n_orbs=len(orbitals)
SO_spin=np.eye(n_spins)
SO_orbital=np.eye(n_orbs)
SO_tensor=-mu*np.kron(SO_spin,SO_orbital)
############ hopping ##############
t={t}
hoppings = [[t, 0, 0, [1,0,0]],
             [t, 0, 0, [0,1,0]],
             [t, 0, 0, [0,0,1]]]
########### impurities #############
V={V}
impurity_loc = {impurity_locations}
impurity_spin = np.arange(n_spins)
impurity_orb = np.arange(n_orbs)
impurities = None #[V, impurity_loc, impurity_spin, impurity_orb]
########### Interaction ###############
U={U}
U_orbs = np.eye(n_orbs)
U_spin = Pauli_x
hubbard_SO = -U*np.kron(U_spin, U_orbs)
#hubbard_tensor =[[hubbard_SO, [0,0,0]],
#                 [hubbard_SO, [1,0,0], True],
#                 [hubbard_SO, [0,1,0], True]]
#n_U=len(hubbard_tensor)
hubbard_U=np.kron(hubbard_SO, np.eye(n_x*n_y*n_z))
############ Mean-field #############
phi, Delta = {phi}, {Delta}

hartree_SO = np.kron(np.diag([phi,phi]),np.eye(n_orbs))
initial_hartree = np.kron(hartree_SO, np.eye(n_x*n_y*n_z))
initial_hartree = np.diag(initial_hartree)

initial_fock = np.diag(np.zeros([n_x*n_y*n_z*n_spins*n_orbs]))

gorkov_SO = np.kron(np.array([[0,+Delta],[-Delta,0]]), np.eye(n_orbs))
initial_gorkov = np.kron(gorkov_SO, np.eye(n_x*n_y*n_z))

external_hartree = None
external_fock = None
external_gorkov= None
############## Energy ###############
epsilon={epsilon}
omega_lim=10*t 
omegas=np.array([0],dtype=complex)+1.0j*epsilon
######## Self-consistent param ########
friction={friction}
max_iterations=100
eps={eps}
T={T}
'''
        with open(confname+'.conf', 'w') as f:
            f.write(temp_conf)
###################################################################
############################# Main ################################
###################################################################
n_spins = 2
n_x=n_y=21
n_list = np.arange(1,21,1)
a = 1.0
mu = 3.57
t = 1.0
V = 0
impurity_locations = [[0,0,0]]
U=2.8
epsilon = 0.1
phi, Delta = 2.08833, 0.488416
eps=0.001
T=0
friction=0

Main()
