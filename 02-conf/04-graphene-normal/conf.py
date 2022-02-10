'''Creates .conf scripts'''
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
def CreateName(n,mu,s,delta,t,Delta,V):
        return f'Normal_{n:d}_{mu:0.2f}_{s:0.2e}_{delta:0.2e}_{t:0.2f}_{Delta:0.2e}_{V:0.2e}'.replace('0.00e+00','0')

def Main():
    ClearOutFolder(out_folder)
    i=0
    for mu in muu:
        confname = os.path.join(out_folder,CreateName(n,mu,s,delta,t,Delta,V))
        temp_conf = f'''from lib import *
####################################################################
########################## Graphene ############################
####################################################################
n_spins={n_spins}
n_x=n_y={n}
a={a}
# 3D
a1=(a/2)*np.array([3,np.sqrt(3),0])
a2=(a/2)*np.array([3,-np.sqrt(3),0])
d3=0.5*np.array([-1,0,0])
zero=np.array([0,0,0])
dimensions=[n_x,n_y,1]
basis=[a1,a2]
######### lattice & sites #########
mu={mu}
orbitals=[zero,d3]
n_orbitals=len(orbitals)
SO_spin=np.eye(n_spins)
SO_orbital=np.eye(n_orbitals)
SO_tensor=-mu*np.kron(SO_spin,SO_orbital)
############ hopping ##############
t={t}
hop_links = [[t, 0, 1, [0,0]],
             [t, 0, 1, [1,0]],
             [t, 0, 1, [0,1]]]
########### impurities #############
V={V}
impurity_loc = {impurity_locations}
impurity_spin = np.arange(n_spins)
impurity_orb = np.arange(n_orbitals)
impurities = [impurity_loc, impurity_spin, impurity_orb]
############# SC state ############

Delta={Delta}
SC_spin=1#1.0j*Pauli_y
SC_orbital=Delta*np.eye(n_orbitals)
SC_tensor=np.kron(SC_spin,SC_orbital)
############## Energy ###############
epsilon={epsilon}
omega_pts={omega_pts}
omega_lim=10*t 
# omegas=np.linspace(-omega_lim,omega_lim,omega_pts,endpoint=True,dtype=complex)+1.0j*epsilon
omegas=np.array([0],dtype=complex)+1.0j*epsilon'''
        with open(confname+'.conf', 'w') as f:
            f.write(temp_conf)
        i+=1
    print(f'{i} .conf files created')
###################################################################
############################# Main ################################
###################################################################
n_spins = 1
n_orbitals = 2
n = 3
a = 1.0
muu = np.linspace(-4.10, 4.10, 1, dtype=float)
s = 0.05
delta = 0.075
t = 1.0
Delta = 0
V =1.21
impurity_locations = [[0,0]]
epsilon = 0.1
omega_pts = 401

Main()
