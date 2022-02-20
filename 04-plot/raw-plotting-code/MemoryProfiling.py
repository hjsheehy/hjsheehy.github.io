from lib_plt import *

import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

logging.debug('importing')
tt=time.perf_counter()

# out_folder='../out/'
out_folder=''

confname='Test'
#
def Conf(
n_spins,
n_orbitals,
n,
mu,
s,
delta,
t,
V,
impurity_locations,
U,
Delta):
    ###################################################################
    ######################### Normal state ############################
    ###################################################################
    onsite_spin=-mu*np.eye(n_spins)
    onsite_orbital=np.eye(n_orbitals)#-s*Pauli_z+delta*Pauli_x
    onsite_tensor=np.kron(onsite_spin,onsite_orbital)    
    
    nn_spin=np.eye(n_spins)
    nn_orbital=-t*np.eye(n_orbitals)
    nn_tensor=np.kron(nn_spin,nn_orbital)

    n_x=n_y=n
    n_sites = n_x * n_y
    dof = n_sites * n_spins * n_orbitals
    ###################################################################
    ########################## Impurities #############################
    ###################################################################
    im_spin=np.eye(n_spins)
    im_orbital=V*np.eye(n_orbitals)
    impurity_tensor=np.kron(im_spin,im_orbital)
    ###################################################################
    ########################## SC state ###############################
    ###################################################################
    SC_spin=1.0j*Pauli_y
    SC_orbital=Delta*np.eye(n_orbitals)    
    SC_tensor=np.kron(SC_spin,SC_orbital)
    ###################################################################
    ############################ Energy ###############################
    ###################################################################
    epsilon=0.1
    omega_pts=401
    omega_lim=10*t 
    omegas=np.linspace(-omega_lim,omega_lim,omega_pts,endpoint=True,dtype=complex)+1.0j*epsilon
    ###################################################################
    ########################### Lattice ###############################
    ###################################################################
    basis_lattice = np.array([[1,0], [0,1]])
    basis_unit_cell = np.array([[0.1, 0.3], [0.3, 0.1], [0.25, 0.1]])
    
    ###################################################################
    ############### Warning: Don't try this at home  ##################
    ###################################################################
    globals().update(locals())
    
##########################################################
########################## Conf ##########################
########################################################## 
n_spins=2
n_orbitals=2
n=35
mu=-0.67
s=0
delta=0
t=1
V=0
impurity_locations=[[0,0]]
U=2.8
Delta=1.2

Conf(n_spins,n_orbitals,n,mu,s,delta,t,V,impurity_locations,U,Delta)
index = Index(dof, n_x, n_y, n_sites,n_spins,n_orbitals)
########################################################################
############################### Main ###################################
########################################################################
# @profile
def Main():
    h0 = H0(onsite_tensor, nn_tensor, impurity_tensor,impurity_locations, n_x, n_y)
    mf = Mean_field(SC_tensor, n_sites)
    ham = H(h0, mf)
    del(h0)
    del(mf)
    gc.collect() # Garbage collector prevents MemoryError
    memory=sys.getsizeof(ham)
    logging.debug('matrix size (MiB):')
    logging.debug(memory/1048576)

    w,v = la.eigh(ham, overwrite_a=True)
    del(ham)
    
    logging.debug('running DOS')
    density_matrix = DM(v, n_x,n_y,index,n_sites,n_spins,n_orbitals)
    dos = DOS(omegas, w, density_matrix)
    del(density_matrix)
    gc.collect()

    density_matrix = ADM(v, n_x,n_y,index,n_sites,n_spins,n_orbitals)
    ados = DOS(omegas, w, density_matrix)
    
    logging.debug('saving')
    with open(out_folder+confname+'.npz', 'wb') as f:
        np.savez(f,
        dos=dos,
        ados=ados)#,
        # executionTime=executionTime,
        # memory=memory)
    return
Main()
# temp=Lattice(n_x, n_y, n_sites, basis_lattice, basis_unit_cell)
# print(np.shape(temp))
logging.debug('Done')
# if __name__ == '__main__':
    # Main()