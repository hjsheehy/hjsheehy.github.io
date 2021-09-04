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
    n_x = int(f.readline().split('#')[0])
    # n_x_bits = f.readline().split('#')[0].split(',')
    # n_x_bits = [int(n_x_bits[0]), int(n_x_bits[1]), int(n_x_bits[2])]
    # n_sites_bits = f.readline().split('#')[0].split(',')
    mu = float(f.readline().split('#')[0])
    s= float(f.readline().split('#')[0])
    delta = float(f.readline().split('#')[0])
    t = float(f.readline().split('#')[0])
    # Delta = float(f.readline().split('#')[0])
    Delta_bits = f.readline().split('#')[0].split(',')
    Delta_bits = [float(Delta_bits[0]), float(Delta_bits[1]), float(Delta_bits[2])]
    V = float(f.readline().split('#')[0])
    impurity_locations = f.readline().split('#')[0]

i=0
def create_name(n_x,n_y,mu,s,delta,t,Delta,V):
        return f'Singlet_varying_Delta_{n_x:d}_{n_y:d}_{mu:0.2f}_{s:0.2e}_{delta:0.2e}_{t:0.2f}_{Delta:0.2e}_{V:0.2e}'.replace('0.00e+00','0')

# for n_x in np.arange(n_sites_bits[0],n_sites_bits[1]+0.01*n_sites_bits[2],n_sites_bits[2],dtype=int):
for Delta in np.linspace(Delta_bits[0],Delta_bits[1],Delta_bits[2],endpoint=True):
    # print(i)
    # i+=1
    n_y=n_x
    confname = os.path.join(fileDir, 'output/'+create_name(n_x,n_y,mu,s,delta,t,Delta,V))
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
onsite_orbital=-mu*np.eye(n_orbitals)#-s*Pauli_z+delta*Pauli_x
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
SC_spin=1.0j*Pauli_y
SC_orbital=Delta*np.eye(n_orbitals)
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