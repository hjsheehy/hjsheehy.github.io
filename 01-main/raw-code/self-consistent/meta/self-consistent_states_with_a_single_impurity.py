'''Creates conf scripts for USC_eigsys'''
import sys
import os
import numpy as np
###########################################################################
fileDir = os.path.dirname(os.path.realpath('__file__'))
print(fileDir)

if len(sys.argv) < 2:
	print('Please supply name of configuration file')
	sys.exit()
    
configuration = sys.argv[1]

with open(configuration,'r') as f:
    n_spins = int(f.readline().split('#')[0])
    n_orbitals = int(f.readline().split('#')[0])
    n = int(f.readline().split('#')[0])
    mu = float(f.readline().split('#')[0])
    s= float(f.readline().split('#')[0])
    delta = float(f.readline().split('#')[0])
    t = float(f.readline().split('#')[0])
    V = float(f.readline().split('#')[0])
    impurity_locations = f.readline().split('#')[0]
    U = float(f.readline().split('#')[0])
    phi_up_plus = float(f.readline().split('#')[0])
    phi_down_plus = float(f.readline().split('#')[0])
    phi_up_minus = float(f.readline().split('#')[0])
    phi_down_minus = float(f.readline().split('#')[0])
    Delta = float(f.readline().split('#')[0])
    varsigma = float(f.readline().split('#')[0])
    friction = float(f.readline().split('#')[0])
    max_iterations = int(f.readline().split('#')[0])
    eps = float(f.readline().split('#')[0])
    T = float(f.readline().split('#')[0])

i=0
def create_name(n,mu,s,delta,t,Delta,V,U,T):
        return f'self-consistent_{n:d}_{mu:0.2f}_{s:0.2e}_{delta:0.2e}_{t:0.2f}_{V:0.2e}_{U:0.2e}_{T:0.2e}'.replace('0.00e+00','0')

impurity_spin_list = ['np.eye(n_spins)', 'np.eye(n_spins)', 'np.array([[1,0],[0,0]])']
impurity_coupling_list = [0, V, V]
impurity_name_list = ['no_impurity', 'impurity', 'mag_impurity']

Pauli_y = np.array([[0,-1.0j],[1.0j,0]])

state_name_list = ['multiorbital_singlet', 'non-unitary_triplet']

for i in range(len(impurity_name_list)):
    for j in range(len(state_name_list)):
        confname = os.path.join(fileDir, 'output/'+state_name_list[j]+'_'+impurity_name_list[i]+'_'+create_name(n,mu,s,delta,t,Delta,V,U,T))
        temp_conf = f'''from lib import *
###################################################################
######################### Normal state ############################
###################################################################
n_spins={n_spins}
n_orbitals={n_orbitals}
n={n}
n_x=n_y=n
mu={mu}
s={s}
delta={delta}
t={t}

onsite_spin=-mu*np.eye(n_spins)
onsite_orbital=np.eye(n_orbitals)-s*Pauli_z+delta*Pauli_x
onsite_tensor=np.kron(onsite_spin,onsite_orbital)    

nn_spin=np.eye(n_spins)
nn_orbital=-t*np.eye(n_orbitals)
nn_tensor=np.kron(nn_spin,nn_orbital)

n_sites = n_x * n_y
dof = n_sites * n_spins * n_orbitals
###################################################################
########################## Impurities #############################
###################################################################
V={impurity_coupling_list[i]}
impurity_locations={impurity_locations}
impurity_spin={impurity_spin_list[i]}
impurity_orbital=V*np.eye(n_orbitals)
impurity_tensor=np.kron(impurity_spin,impurity_orbital)
###################################################################
######################### Interaction #############################
###################################################################
U={U}
U_orbital = np.eye(n_orbitals)
U_spin = Pauli_x
Hubbard_tensor = -U*np.kron(np.eye(n_sites), np.kron( U_spin, U_orbital))
###################################################################
######################### Mean-field ##############################
###################################################################
Fock = np.zeros([n_x,n_y,n_spins,n_spins,n_orbitals,n_orbitals])

Hartree = np.zeros([n_x,n_y,n_spins,n_orbitals])
Hartree[:,:,0,0]={phi_up_plus}*np.ones([n_x,n_y])
Hartree[:,:,1,0]={phi_down_plus}*np.ones([n_x,n_y])
Hartree[:,:,0,1]={phi_up_minus}*np.ones([n_x,n_y])
Hartree[:,:,1,1]={phi_down_minus}*np.ones([n_x,n_y])

Delta={Delta}
varsigma={varsigma}
sym = np.array([[1+varsigma, 0], [0, 1-varsigma]])
anti_sym = 1.0j*Pauli_y
Gorkov_spin_orbital = {['np.kron(anti_sym, sym)', 'np.kron(sym, anti_sym)'][j]}
Gorkov_tensor=Delta*np.kron(np.eye(n_sites), Gorkov_spin_orbital)

friction={friction}
max_iterations={max_iterations}
eps={eps}
T={T}
###################################################################
############################ Energy ###############################
###################################################################
Include_DOS=True
epsilon=0.1
omega_pts=401
omega_lim=10*t 
omegas=np.linspace(-omega_lim,omega_lim,omega_pts,endpoint=True,dtype=complex)+1.0j*epsilon'''
        with open(confname+'.conf', 'w') as f:
            f.write(temp_conf)
print('Done')