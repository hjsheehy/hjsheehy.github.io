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
def CreateName(n,mu,t,V,U):
    return f'Normal_{n:d}_{mu:0.2f}_{t:0.2f}_{V:0.2e}_{U:0.2e}'.replace('0.00e+00','0')

def Main():
    ClearOutFolder(out_folder)
    i=0
    for mu in muu:
        confname = os.path.join(out_folder,CreateName(n,mu,t,V,U))
        temp_conf = f'''from lib import *
####################################################################
######################## Self-consistent ###########################
####################################################################
n_spins={n_spins}
n_x=n_y={n}
n_ph={n_ph}
a={a}
a1=(a/2)*np.array([3,np.sqrt(3),0])
a2=(a/2)*np.array([3,-np.sqrt(3),0])
d3=0.5*np.array([-1,0,0])
zero=np.array([0,0,0])
dimensions=[n_x,n_y,1]
basis=[a1,a2]
######### lattice & sites #########
mu={mu}
orbitals=[zero]
n_orbs=len(orbitals)
SO_spin=np.eye(n_spins)
SO_orbs=np.eye(n_orbs)
SO_tensor=-mu*np.kron(SO_spin,SO_orbs)
############ hopping ##############
t={t}
hop_links = [[t, 0, 0, [1,0]],
             [t, 0, 0, [0,1]]]
########### impurities #############
V={V}
impurity_loc = {impurity_locations}
impurity_spin = np.arange(n_spins)
impurity_orb = np.arange(n_orbs)
impurities = [impurity_loc, impurity_spin, impurity_orb]
########### Interaction ###############
U={U}
U_orbs = np.eye(n_orbs)
U_spin = Pauli_x
Hubbard_tensor = -U*np.kron(np.eye(np.prod(dimensions)), np.kron( U_spin, U_orbs))
############ Mean-field #############
phiUp, phiDown, Delta= {phiUp}, {phiDown}, {Delta}

Fock = np.zeros([dimensions[0],dimensions[1],dimensions[2],n_spins,n_spins,n_orbs,n_orbs])

Hartree = np.zeros([dimensions[0],dimensions[1],dimensions[2],n_spins,n_orbs])
Hartree[:,:,:,0,:]=-phiUp*np.ones([dimensions[0],dimensions[1],dimensions[2],n_orbs])
Hartree[:,:,:,1,:]=-phiDown*np.ones([dimensions[0],dimensions[1],dimensions[2],n_orbs])

Gorkov= np.zeros([dimensions[0],dimensions[1],dimensions[2],n_spins,n_spins,n_orbs,n_orbs])
Gorkov[:,:,:,0,1,:,:]=+Delta*np.ones([dimensions[0],dimensions[1],dimensions[2],n_orbs,n_orbs])
Gorkov[:,:,:,1,0,:,:]=-Delta*np.ones([dimensions[0],dimensions[1],dimensions[2],n_orbs,n_orbs])
############## Energy ###############
epsilon={epsilon}
omega_pts={omega_pts}
omega_lim=10*t 
omegas=np.linspace(-omega_lim,omega_lim,omega_pts,endpoint=True,dtype=complex)+1.0j*epsilon
######## Self-consistent param ########
friction={friction}
max_iterations={max_iterations}
eps={eps}
T={T}'''
        with open(confname+'.conf', 'w') as f:
            f.write(temp_conf)
        i+=1
    print(f'{i} .conf files created')
###################################################################
############################# Main ################################
###################################################################
n_spins = 2
n = 31
n_ph=2
a = 1.0
#muu = np.linspace(-4.10, 4.10, 3, dtype=float)
muu=[-3.51]
t = 1.0
V =1.21
U=2.23
phiUp=1.3
phiDown=1.7
Delta=0.9
impurity_locations = [[0,0]]

epsilon = 0.1
omega_pts = 401

friction=0.1
max_iterations=20
eps=0.1
T=0

Main()
