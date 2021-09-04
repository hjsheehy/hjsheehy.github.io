'''Creates conf scripts for USC_eigsys'''
import sys
import os
import numpy as np
import random
###########################################################################
fileDir = os.path.dirname(os.path.realpath('__file__'))

if len(sys.argv) < 2:
	print('Please supply name of configuration file')
	sys.exit()
    
configuration = sys.argv[1]

with open(configuration,'r') as f:
    n_spins = int(f.readline().split('#')[0])
    n_orbitals = int(f.readline().split('#')[0])
    n_x = int(f.readline().split('#')[0])
    # n_x_bits = f.readline().split('#')[0].split(',')
    # n_x_bits = [int(n_x_bits[0]), int(n_x_bits[1]), int(n_x_bits[2])]
    mu= float(f.readline().split('#')[0])
    # mu_bits = f.readline().split('#')[0].split(',')
    # mu_bits = [float(mu_bits[0]), float(mu_bits[1]), float(mu_bits[2])]
    s= float(f.readline().split('#')[0])
    delta = float(f.readline().split('#')[0])
    t = float(f.readline().split('#')[0])
    Delta = float(f.readline().split('#')[0])
    eta_bits = f.readline().split('#')[0].split(',')
    eta = [complex(eta_bits[0]), complex(eta_bits[1]), complex(eta_bits[2])]
    V = float(f.readline().split('#')[0])
    # V_bits = f.readline().split('#')[0].split(',')
    # V_bits = [float(V_bits[0]), float(V_bits[1]), float(V_bits[2])]
    n_imp = int(f.readline().split('#')[0])
    n_trials = int(f.readline().split('#')[0])

i=0
def create_name(n_x,n_y,mu,s,delta,t,Delta,eta,V,trial):
        return f'Random_impurities_triplet_{trial:d}_{n_x:d}_{n_y:d}_{mu:0.2f}_{s:0.2e}_{delta:0.2e}_{t:0.2f}_{Delta:0.2e}_{eta[0]}_{eta[1]}_{eta[2]}_{V:0.2e}'.replace('0.00e+00','0')

for trial in range(n_trials):
    n_y=n_x
    confname = os.path.join(fileDir, 'output/'+create_name(n_x,n_y,mu,s,delta,t,Delta,eta,V,trial))
    impurity_locations=[[random.randint(0,n_x),random.randint(0,n_y)] for i in range(n_imp)]
    temp_conf = f'''from lib import *
###################################################################
######################### Normal state ############################
###################################################################
n_spins={n_spins}
n_orbitals={n_orbitals}
n_x={n_x}
n_y={n_y}
mu={mu}
s={s}
delta={delta}
t={t}
onsite_spin=np.eye(n_spins)
onsite_orbital=-mu*np.eye(n_orbitals)-s*Pauli_z+delta*Pauli_x
nn_spin=np.eye(n_spins)
nn_orbital=-t*np.eye(n_orbitals)
n_sites = n_x * n_y
dof = n_sites * n_spins * n_orbitals
###################################################################
########################## Impurities #############################
###################################################################
V={V}
impurity_locations={impurity_locations}
im_spin=np.eye(n_spins)
im_orbital=V*np.eye(n_orbitals)
###################################################################
########################### SC state ##############################
###################################################################
Delta={Delta}
eta=np.array({eta})
d=Delta*eta
SC_spin=1.0j*np.matmul((np.einsum('k,ijk->ij',d,Pauli_vec)),Pauli_y)
SC_orbital=1.0j*Pauli_y
###################################################################
############################ Energy ###############################
###################################################################
epsilon=0.1
omega_pts=401
omega_lim=10*t 
omegas=np.linspace(-omega_lim,omega_lim,omega_pts,endpoint=True,dtype=complex)+1.0j*epsilon'''
    with open(confname+'.conf', 'w') as f:
        f.write(temp_conf)
print('Done')