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
    U_bits = f.readline().split('#')[0].split(',')
    U_bits = [float(U_bits[0]), float(U_bits[1]), float(U_bits[2])]
    phiUp = float(f.readline().split('#')[0])
    phiDown = float(f.readline().split('#')[0])
    Delta = float(f.readline().split('#')[0])
    friction = float(f.readline().split('#')[0])
    max_iterations = int(f.readline().split('#')[0])
    eps = float(f.readline().split('#')[0])
    T_bits = f.readline().split('#')[0].split(',')
    T_bits = [float(T_bits[0]), float(T_bits[1]), float(T_bits[2])]

i=0
def create_name(n,mu,s,delta,t,Delta,V,U,T):
        return f'Self-consistent_singlet_{n:d}_{mu:0.2f}_{s:0.2e}_{delta:0.2e}_{t:0.2f}_{V:0.2e}_{U:0.2e}_{T:0.2e}'.replace('0.00e+00','0')

for U in np.linspace(U_bits[0],U_bits[1],U_bits[2],endpoint=True):
    for T in np.linspace(T_bits[0],T_bits[1],T_bits[2],endpoint=True):
        confname = os.path.join(fileDir, 'output/'+create_name(n,mu,s,delta,t,Delta,V,U,T))
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
onsite_orbital=np.eye(n_orbitals)#-s*Pauli_z+delta*Pauli_x
onsite_tensor=np.kron(onsite_spin,onsite_orbital)    

nn_spin=np.eye(n_spins)
nn_orbital=-t*np.eye(n_orbitals)
nn_tensor=np.kron(nn_spin,nn_orbital)

n_sites = n_x * n_y
dof = n_sites * n_spins * n_orbitals
###################################################################
########################## Impurities #############################
###################################################################
V={V}
impurity_locations={impurity_locations}
im_spin=np.eye(n_spins)
im_orbital=V*np.eye(n_orbitals)
impurity_tensor=np.kron(im_spin,im_orbital)
###################################################################
######################### Interaction #############################
###################################################################
U={U}
U_orbital = np.eye(n_orbitals)
U_spin = np.eye(n_spins) #Pauli_x
Hubbard_tensor = -U*np.kron(np.eye(n_sites), np.kron( U_spin, U_orbital))
###################################################################
######################### Mean-field ##############################
###################################################################
phiUp, phiDown, Delta= {phiUp}, {phiDown}, {Delta}

Fock = np.zeros([n_x,n_y,n_spins,n_spins,n_orbitals,n_orbitals])

Hartree = np.zeros([n_x,n_y,n_spins,n_orbitals])
Hartree[:,:,0,:]=-phiUp*np.ones([n_x,n_y,n_orbitals])
Hartree[:,:,1,:]=-phiDown*np.ones([n_x,n_y,n_orbitals])

Gorkov= np.zeros([n_x,n_y,n_spins,n_spins,n_orbitals,n_orbitals])
Gorkov[:,:,0,1,:,:]=+Delta*np.ones([n_x,n_y,n_orbitals,n_orbitals])
Gorkov[:,:,1,0,:,:]=-Delta*np.ones([n_x,n_y,n_orbitals,n_orbitals])

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