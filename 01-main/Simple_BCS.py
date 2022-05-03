from lib import *

def main():
    A=Atom([0,0],'A')
    A.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    bdg=BogoliubovdeGennes(lattice_vectors,'Simple-BCS')
    bdg.add_atom(A)
    bdg.n_spins=2
 
    if alpha:
        bdg.set_kpts([n_cells,n_cells])
    else:
        bdg.cut(n_cells, axes=0, glue_edgs=True)
        bdg.cut(n_cells, axes=1, glue_edgs=True)

    bdg.set_onsite(-mu)
    bdg.set_hopping(-t,hop_vector=[1,0],label='$t$')
    bdg.set_hopping(-t,hop_vector=[0,1],label='$t$')

    # bdg.add_impurities(V,[0,0])
    bdg.set_hartree(rho)
    bdg.set_fock(phi,spin_i='up',spin_f='dn')
    bdg.set_gorkov(chi*1j*Pauli_y)
    
    bdg.set_hubbard_u(-U,spin_i='up',spin_f='dn')

    # bdg._set_mean_field_hamiltonian()
    # hamiltonian = np.real(bdg._hamiltonian)
    # plt.imshow(hamiltonian)
    # plt.show()
    # exit()

    # _print=True
    # bdg.record_hartree(location=[0,0], spin='up', _print=_print)
    # bdg.record_hartree(location=[0,0], spin='down', _print=_print)
    # bdg.record_fock(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='down', _print=_print)
    # bdg.record_fock(location_i=[0,0], location_f=[0,0], spin_i='down', spin_f='up', _print=_print)
    # bdg.record_gorkov(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='down',_print=_print)
    # bdg.record_gorkov(location_i=[0,0], location_f=[0,0], spin_i='down', spin_f='up',_print=_print)

    bdg.self_consistent_calculation(friction=0.2, max_iterations=400, absolute_convergence_factor=0.0001)
    # bdg.solve()
    
    print(bdg._gorkov)
    print('#######')

    energy_interval=np.linspace(-4,4,601)
    resolution=0.02

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
mu=3.2
t=1
U=3.6
rho=5.2
phi=0.1
chi=2.2
# rho=0
# phi=0
# chi=0
# V=1.21
n_cells=11

for alpha in [True, False]:
    # greens_function_xy, greens_function_xq, greens_function_kq, bdg = main()
    greens_function_kq = main()

    # k-space
    # fig,ax = plt.subplots(1,1)
    # greens_function_kq.plot_ldos(ax, energy=0, anomalous=False)
# plt.show()
