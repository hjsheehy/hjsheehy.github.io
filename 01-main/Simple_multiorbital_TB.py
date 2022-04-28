from lib import *

def main(alpha):
    A=Atom([0,0],'A')
    B=Atom([0.5,0],'B')
    A.add_orbital('s')
    B.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    tb=TightBinding(lattice_vectors,'Simple multiorbital TB')
    tb.add_atom(A)
    tb.add_atom(B)
    tb.n_spins=1
    
    if alpha:
        tb.set_kpts([n_cells,n_cells])
    else:
        tb.cut(n_cells, axes=0, glue_edgs=True)
        tb.cut(n_cells, axes=1, glue_edgs=True)

    tb.set_onsite(-mu+s,atom='A')
    tb.set_onsite(-mu-s,atom='B')

    tb.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B',label='$v$')
    tb.set_hopping(-w,hop_vector=[-1,0],atom_i='A',atom_f='B',label='$w$')
    tb.set_hopping(-td,hop_vector=[0,1],atom_i='B',atom_f='A',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[0,-1],atom_i='B',atom_f='A',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[1,1],atom_i='B',atom_f='A',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[1,-1],atom_i='B',atom_f='A',label='$t_d$')
    
    # impurity_wall = [[0,i] for i in range(n_cells)]
    # tb.add_impurities(V,impurity_wall)

    # tb.add_impurities(V,[0,0])

    # tb.set_hartree(rho)
    # tb.set_fock(phi,atom_i='A',atom_f='B')
    # tb.set_fock(phi,atom_i='B',atom_f='A',hop_vector=[1,0])
    # tb.set_gorkov(chi,atom_i='A',atom_f='B')
    # tb.set_gorkov(chi,atom_i='B',atom_f='A',hop_vector=[1,0])

    # tb.set_hubbard_u(-Uv,atom_i='A',atom_f='B',hop_vector=[0,0],orbital_i='s',orbital_f='s')
    # tb.set_hubbard_u(Uw,atom_i='B',atom_f='A',hop_vector=[1,0])

    # tb.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B', orbital_i='s', orbital_f='s', spin_i=0, spin_f=0,_print=False)
    # tb.record_gorkov(location_i=[0,0], location_f=[1,0], atom_i='B', atom_f='A', orbital_i='s', orbital_f='s', spin_i=0, spin_f=0,_print=False)

    # tb.self_consistent_calculation(friction=0., max_iterations=2, absolute_convergence_factor=0.00001)
    tb.solve()
    
    energy_interval=np.linspace(-4,4,601)
    resolution=0.05

    # greens_function_xy=GreensFunction(tb,energy_interval,resolution, k_axes=None)

    # greens_function_xq=GreensFunction(tb,energy_interval,resolution, k_axes=[1])
<<<<<<< HEAD

    # print(np.real(tb._hamiltonian))
    # print('####')
=======
    # print(tb.kpts)
    # k=0
    # tmp=-v-2*w*np.cos(k)-4*2*td*np.cos(k)
    # print(tmp)
    # print(np.real(tb._hamiltonian))
    # print('####')
    # exit()
>>>>>>> e3f6c4460c17734a30563a387a3ab2800fa31082
    greens_function_kq=GreensFunction(tb,energy_interval,resolution, k_axes=[0,1])

    del tb.eigenvectors
    del tb.eigenvalues

    return greens_function_kq
    # return greens_function_xy, greens_function_xq, greens_function_kq

# Plotting:

#############################################################################
################################# Main ######################################
#############################################################################
<<<<<<< HEAD
mu=-3.2
s=0
w=1.2
v=0.6
td=0#.9
V=0
n_cells=11
=======
mu=-3.8
s=0
w=1.2
v=0.6
# w=v=0
td=0.9
# td=0
V=0
n_cells=21
>>>>>>> e3f6c4460c17734a30563a387a3ab2800fa31082

fig,ax = plt.subplots(2,2)
greens_function_kq=[]
<<<<<<< HEAD
=======
title=['Continuum model','Lattice model']
>>>>>>> e3f6c4460c17734a30563a387a3ab2800fa31082
for i,alpha in enumerate([True,False]):
    # greens_function_xy, greens_function_xq, greens_function_kq = main()
    greens_function_kq.append(main(alpha))

    #real space
    # fig,ax = plt.subplots(1,1)
    # greens_function_xy.plot_ldos(ax, energy=0)
    # plt.show()

    # k-space
    greens_function_kq[i].plot_ldos(ax[0,i], energy=0,atom='A')
    greens_function_kq[i].plot_ldos(ax[1,i], energy=0,atom='B')
    ax[0,i].set_title(title[i])
    ax[1,i].set_title('')
plt.show()
