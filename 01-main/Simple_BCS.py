from lib import *

def main():
    A=Atom([0,0],'A')
    A.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    bdg=BogoliubovdeGennes(lattice_vectors,'Simple-BCS')
    bdg.add_atom(A)
    bdg.n_spins=2
 
    if bulk_calculation:
        bdg.set_kpts([n_cells,n_cells])
    else:
        bdg.cut(n_cells, axes=0, glue_edgs=True)
        bdg.cut(n_cells, axes=1, glue_edgs=True)

    bdg.set_onsite(-mu)
    bdg.set_hopping(-t,hop_vector=[1,0],label='$t$')
    bdg.set_hopping(-t,hop_vector=[0,1],label='$t$')

    # bdg.add_impurities(V,[0,0])
    bdg.set_hartree(rho_up,spin='up')
    bdg.set_hartree(rho_dn,spin='dn')
    bdg.set_fock(phi,spin_i='up',spin_f='dn')
    bdg.set_gorkov(chi*1j*Pauli_y)
    
    # bdg.set_hubbard_u(U,spin_i='up',spin_f='dn')
    bdg.set_hubbard_u(U*Pauli_x)

    # bdg._set_mean_field_hamiltonian()
    # hamiltonian = np.real(bdg._hamiltonian)
    # plt.imshow(hamiltonian)
    # plt.show()
    # exit()

    bdg.record_hartree(location=[0,0], spin='up', _print=_print)
    bdg.record_hartree(location=[0,0], spin='down', _print=_print)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='down', _print=_print)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], spin_i='down', spin_f='up', _print=_print)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='down',_print=_print)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], spin_i='down', spin_f='up',_print=_print)
        
    if print_free_energy:
        bdg.print_V=True
        bdg.print_V_mf=True
        bdg.print_Eg=True
        bdg.print_free_energy=True

    bdg.self_consistent_calculation(friction=friction, max_iterations=max_iterations, absolute_convergence_factor=0.0001)
    # bdg.solve()
    
    energy_interval=np.linspace(-4,4,601)
    resolution=0.02

    # greens_function_xy=GreensFunction(bdg,energy_interval,resolution, k_axes=None)

    # greens_function_xq=GreensFunction(bdg,energy_interval,resolution, k_axes=[1])

    # greens_function_kq=GreensFunction(bdg,energy_interval,resolution, k_axes=[0,1])

    del bdg.eigenvectors
    del bdg.eigenvalues

    # with open(DATA+'.npz', 'wb') as f:
    #     cPickle.dump([greens_function_xy, greens_function_xq, greens_function_kq, bdg], f)
    return bdg
    # return greens_function_kq
    # return greens_function_xy, greens_function_xq, greens_function_kq, bdg

def plot_iterations(bdg):

    markers=['o','+','^','x','.']
    s=3

    # hartree_A=bdg.hartree(atom='A')[0,0]
    # hartree_B=bdg.hartree(atom='B')[0,0]
    # fock_v=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])[0,0]
    # gorkov_v=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])[0,0]
    # free_energy=bdg.free_energy
    
    # exit()
    # gorkov_w=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[-1,0])

    fig, ax1  = plt.subplots(1,1,sharex='col')

    # color = 'tab:black'
    ax1.tick_params(axis='y')
    ax1.plot(bdg._fock_iterations[0],c='k',marker=markers[2],markersize=s,label=f'$\Phi$')
    ax1.plot(bdg._gorkov_iterations[0],c='cyan',marker=markers[3],markersize=s,label=f'$\Delta$')
    ax1.plot(bdg._hartree_iterations[0],c='b',marker=markers[0],markersize=s,label=r'$\phi_\uparrow$')
    ax1.plot(bdg._hartree_iterations[1],c='g',marker=markers[1],markersize=s,label=r'$\phi_\downarrow$')
    ax1.legend()

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('Free energy', color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.plot(bdg.free_energy,c='r',marker=markers[4],markersize=s,label=f'Free energy')
    # ax2.plot(bdg.Eg,c='b',marker=markers[4],markersize=s,label=f'Eg')
    # ax2.plot(bdg.V,c='g',marker=markers[4],markersize=s,label=f'V')
    # ax2.plot(bdg.V_mf,c='c',marker=markers[4],markersize=s,label=r'V_{mf}')

    # plt.plot(np.real(bdg._gorkov_iterations[1]))
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    xlabel=r'Iterations'
    ylabel=r'Amplitude of fields'
    title=r'BCS theory'
    fig.suptitle(title)
    fig.supxlabel(xlabel)
    fig.supylabel(ylabel)
    plt.tight_layout()
    
#############################################################################
################################# Main ######################################
#############################################################################
mu=-2.7
t=1
U=2.8
rho_up=3.4
rho_dn=1.4
phi=0
chi=1.0
n_cells=21
max_iterations=100
friction=0.4

_print=True
print_free_energy=True

# for alpha in [True, False]:
    # greens_function_xy, greens_function_xq, greens_function_kq, bdg = main()
    # greens_function_kq = main()

    # k-space
    # fig,ax = plt.subplots(1,1)
    # greens_function_kq.plot_ldos(ax, energy=0, anomalous=False)

    # iterations
    # fig,ax = plt.subplots(1,1)
# plt.show()

bulk_calculation = False
bdg = main()
# exit()
plot_iterations(bdg)
plt.show()
