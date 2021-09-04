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
    # n_bits = f.readline().split('#')[0].split(',')
    # n_bits = [int(n_bits[0]), int(n_bits[1]), int(n_bits[2])]
    # mu= float(f.readline().split('#')[0])
    mu_bits = f.readline().split('#')[0].split(',')
    mu_bits = [float(mu_bits[0]), float(mu_bits[1]), int(mu_bits[2])]
    s= float(f.readline().split('#')[0])
    delta = float(f.readline().split('#')[0])
    t = float(f.readline().split('#')[0])
    Delta = float(f.readline().split('#')[0])
    V = float(f.readline().split('#')[0])
    # V_bits = f.readline().split('#')[0].split(',')
    # V_bits = [float(V_bits[0]), float(V_bits[1]), float(V_bits[2])]
    impurity_locations = f.readline().split('#')[0]

i=0
def create_name(n,mu,s,delta,t,Delta,V):
        return f'Normal_{n:d}_{mu:0.2f}_{s:0.2e}_{delta:0.2e}_{t:0.2f}_{Delta:0.2e}_{V:0.2e}'.replace('0.00e+00','0')

# for n in np.arange(n_bits[0],n_bits[1]+0.01*n_bits[2],n_bits[2],dtype=int):
for mu in np.linspace(mu_bits[0],mu_bits[1],mu_bits[2],dtype=float):
# for V in np.linspace(V_bits[0],V_bits[1],V_bits[2],dtype=float):
    confname = os.path.join(fileDir, 'out/'+create_name(n,mu,s,delta,t,Delta,V))
    temp_conf = f'''from lib import *
###################################################################
######################### Normal state ############################
###################################################################
n_spins={n_spins}
n_orbitals={n_orbitals}
n={n}
mu={mu}
s={s}
delta={delta}
t={t}
onsite_spin=np.eye(n_spins)
onsite_orbital=-mu*np.eye(n_orbitals)-s*Pauli_z+delta*Pauli_x
onsite_tensor=np.kron(onsite_spin,onsite_orbital)    

nn_spin=np.eye(n_spins)
nn_orbital=-t*np.eye(n_orbitals)
nn_tensor_x=np.kron(nn_spin,nn_orbital)
nn_tensor_y=nn_tensor_x


n_x=n_y=n
n_sites = n_x * n_y
dof = n_sites * n_spins * n_orbitals
###################################################################
########################## Impurities #############################
###################################################################
V={V}
impurity_locations={impurity_locations}
impurity_spin=np.eye(n_spins)
impurity_orbital=V*np.eye(n_orbitals)
impurity_tensor=np.kron(impurity_spin,impurity_orbital)
###################################################################
########################### SC state ##############################
###################################################################
Delta={Delta}
SC_spin=1#1.0j*Pauli_y
SC_orbital=Delta*np.eye(n_orbitals)
SC_tensor=np.kron(SC_spin,SC_orbital)
###################################################################
############################ Energy ###############################
###################################################################
epsilon=0.1
omega_pts=401
omega_lim=10*t 
# omegas=np.linspace(-omega_lim,omega_lim,omega_pts,endpoint=True,dtype=complex)+1.0j*epsilon
omegas=np.array([0],dtype=complex)+1.0j*epsilon'''
    with open(confname+'.conf', 'w') as f:
        f.write(temp_conf)
print('Done')
