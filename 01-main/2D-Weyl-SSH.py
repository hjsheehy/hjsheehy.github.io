from lib import *

filename=sys.argv[0].split('.')[0]
data=DATA+filename+'.npz'

def main():
    A=Atom([0,0],'A')
    B=Atom([0.5,0],'B')
    A.add_orbital('s')
    B.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    bdg=BogoliubovdeGennes(lattice_vectors,'2D-Weyl-SSH')
    bdg.add_atom(A)
    bdg.add_atom(B)
    bdg.n_spins=1
    mu=0.0
    s=0.0
    td=0.9
    v=0.6
    w=1.2
    Uv=1
    # Uw=1
    rho=0
    phi=0
    chi=0.
    V=1.21
    n_cells=41
    bdg.cut(n_cells, axes=0, glue_edgs=False)
    bdg.cut(n_cells, axes=1, glue_edgs=True)
    bdg.set_onsite(-mu+s,atom='A')
    bdg.set_onsite(-mu-s,atom='B')

    def papers_convention():
        bdg.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B',label='$v$')
        bdg.set_hopping(-w,hop_vector=[-1,0],atom_i='A',atom_f='B',label='$w$')
        bdg.set_hopping(-td,hop_vector=[0,1],atom_i='B',atom_f='A',label='$t_d$')
        bdg.set_hopping(-td,hop_vector=[0,-1],atom_i='B',atom_f='A',label='$t_d$')
        bdg.set_hopping(-td,hop_vector=[1,1],atom_i='B',atom_f='A',label='$t_d$')
        bdg.set_hopping(-td,hop_vector=[1,-1],atom_i='B',atom_f='A',label='$t_d$')

    def my_convention():
        bdg.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B',label='$v$')
        bdg.set_hopping(-w,hop_vector=[1,0],atom_i='B',atom_f='A',label='$w$')
        bdg.set_hopping(-td,hop_vector=[0,1],atom_i='A',atom_f='B',label='$t_d$')
        bdg.set_hopping(-td,hop_vector=[0,1],atom_i='B',atom_f='A',label='$t_d$')
        bdg.set_hopping(-td,hop_vector=[-1,1],atom_i='A',atom_f='B',label='$t_d$')
        bdg.set_hopping(-td,hop_vector=[1,1],atom_i='B',atom_f='A',label='$t_d$')

    papers_convention()
    # my_convention()
    
    # impurity_wall = [[0,i] for i in range(n_cells)]
    # bdg.add_impurities(V,impurity_wall)

    # bdg.add_impurities(V,[0,4])

    # bdg.set_hartree(rho)
    # bdg.set_fock(phi,atom_i='A',atom_f='B')
    # bdg.set_fock(phi,atom_i='B',atom_f='A',hop_vector=[1,0])
    bdg.set_gorkov(chi,atom_i='A',atom_f='B')
    # bdg.set_gorkov(chi,atom_i='B',atom_f='A',hop_vector=[1,0])

    bdg.set_hubbard_u(-Uv,atom_i='A',atom_f='B',hop_vector=[0,0],orbital_i='s',orbital_f='s')
    # bdg.set_hubbard_u(Uw,atom_i='B',atom_f='A',hop_vector=[1,0])

    # bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B', orbital_i='s', orbital_f='s', spin_i=0, spin_f=0,_print=False)
    # bdg.record_gorkov(location_i=[0,0], location_f=[1,0], atom_i='B', atom_f='A', orbital_i='s', orbital_f='s', spin_i=0, spin_f=0,_print=False)

    # bdg.self_consistent_calculation(friction=0., max_iterations=2, absolute_convergence_factor=0.00001)
    bdg.solve()

    energy_interval=np.linspace(-4,4,601)
    resolution=0.05

    greens_function_xy=GreensFunction(bdg,energy_interval,resolution, k_axes=None)

    greens_function_xq=GreensFunction(bdg,energy_interval,resolution, k_axes=[1])

    greens_function_kq=GreensFunction(bdg,energy_interval,resolution, k_axes=[0,1])

    with open(data, 'wb') as f:
        cPickle.dump([greens_function_xy, greens_function_xq, greens_function_kq, bdg], f)
    return greens_function_xy, greens_function_xq, greens_function_kq, bdg

# Plotting:

def iterations(bdg):
    gorkov_v=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])
    # gorkov_w=bdg.gorkov(atom_i='B', atom_f='A', hop_vector=[1,0])

    plt.plot(np.real(bdg._gorkov_iterations[0]))
    plt.plot(np.real(bdg._gorkov_iterations[1]))
    plt.show()
    plt.close()
    
def unit_cell(model):
    fig, ax = plt.subplots(1, 1)
    model.plot_unit_cell(fig, ax, atoms='all', s=100)
    plt.show()
    plt.close()

def ldos_each_atom(greens_function_xy):
    ldos=greens_function_xy.local_density_of_states(energy=0, atom='A')
    ldos=np.fft.fftshift(ldos.T)
    plt.imshow(ldos)
    plt.show()
    plt.close()

    ldos=greens_function_xy.local_density_of_states(energy=0, atom='B')
    ldos=np.fft.fftshift(ldos.T)
    plt.imshow(ldos)
    plt.show()
    plt.close()

def real_space(greens_function):

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, axes=['resolved','integrated'],omega_min=-2,omega_max=2,vmin=0,vmax=10)

    plt.show()
    plt.close()

def k_space(greens_function):

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, axes=['integrated','resolved'],omega_min=-2,omega_max=2,vmin=0,vmax=80)

    plt.show()
    plt.close()

def majorana_fermi_arc(greens_function):

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, axes=['integrated','resolved'], omega_min=-2,omega_max=2,vmin='default',vmax='default')

    plt.show()
    plt.close()

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, energy=0, axes=['resolved',11], omega_min=0,omega_max='default',vmin=0,vmax=6,label='Majorana')
    ax = greens_function.plot_spectrum(ax, energy=0.5, axes=['resolved',11], omega_min=0,omega_max='default',vmin=0,vmax=6,label='Bogoloiubov-Fermi arc')
    ax.set_title('')
    ax.legend()

    plt.show()
    plt.close()
#############################################################################
################################# Main ######################################
#############################################################################
# greens_function_xy, greens_function_xq, greens_function_kq, bdg = main()

[greens_function_xy, greens_function_xq, greens_function_kq, bdg] = np.load(data, allow_pickle=True)

# plot_itartions(bdg)
# unit_cell(bdg)
# ldos_each_atom(greens_function_xy)
# real_space(greens_function_xy)
# k_space(greens_function_kq)
majorana_fermi_arc(greens_function_xq)
