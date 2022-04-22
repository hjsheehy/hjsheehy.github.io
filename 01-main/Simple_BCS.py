from lib import *

def main():
    A=Atom([0,0],'A')
    A.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    bdg=BogoliubovdeGennes(lattice_vectors,'Simple-BCS')
    bdg.add_atom(A)
    bdg.n_spins=2
    
    bdg.set_kpts([n_cells,n_cells])

    # bdg.cut(n_cells, axes=0, glue_edgs=True)
    # bdg.cut(n_cells, axes=1, glue_edgs=True)
    bdg.set_onsite(-mu)
    bdg.set_hopping(-t,hop_vector=[1,0],label='$t$')
    bdg.set_hopping(-t,hop_vector=[0,1],label='$t$')

    # bdg.add_impurities(V,[0,0])

    bdg.set_hartree(rho)
    bdg.set_fock(phi,spin_i='up',spin_f='dn')
    bdg.set_gorkov(chi,spin_i='up',spin_f='dn')

    bdg.set_hubbard_u(-U,spin_i='up',spin_f='dn')

    # bdg.record_hartree(location=[0,0], atom='A', _print=False)
    # bdg.record_hartree(location=[0,0], atom='B', _print=False)
    # bdg.record_fock(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B',_print=False)
    # bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B',_print=False)
    # bdg.record_gorkov(location_i=[0,0], location_f=[1,0], atom_i='B', atom_f='A', orbital_i='s', orbital_f='s', spin_i=0, spin_f=0,_print=False)

    bdg.self_consistent_calculation(friction=0.2, max_iterations=400, absolute_convergence_factor=0.00001)

    energy_interval=np.linspace(-4,4,601)
    resolution=0.05

    # greens_function_xy=GreensFunction(bdg,energy_interval,resolution, k_axes=None)

    # greens_function_xq=GreensFunction(bdg,energy_interval,resolution, k_axes=[1])

    greens_function_kq=GreensFunction(bdg,energy_interval,resolution, k_axes=[0,1])

    del bdg.eigenvectors
    del bdg.eigenvalues

    # with open(DATA+'.npz', 'wb') as f:
    #     cPickle.dump([greens_function_xy, greens_function_xq, greens_function_kq, bdg], f)
    return greens_function_kq
    # return greens_function_xy, greens_function_xq, greens_function_kq, bdg

#############################################################################
################################# Main ######################################
#############################################################################
mu=-3.7
t=1
U=3.6
rho=-2.2
phi=0.1
chi=1.2
V=0
n_cells=11

main()
