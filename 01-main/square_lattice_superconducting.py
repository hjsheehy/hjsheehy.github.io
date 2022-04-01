from lib import *
A=Atom([0,0],'A')
A.add_orbital('s')
lattice_vectors=[[1,0],[0,1]]
bdg=BogoliubovdeGennes(lattice_vectors,'Sr2RuO4')
bdg.add_atom(A)
bdg.n_spins=2
mu=-3.72
t=1
U=2.18
rho=0#2.21
phi=0#0.32
chi=0#1.42
V=2.21
n_cells=11
bdg.cut_piece(n_cells, [0,1])
bdg.set_onsite(-mu,orbital='s')
bdg.set_hopping(-t,hop_vector=[1,0])
bdg.set_hopping(-t,hop_vector=[0,1])
bdg.add_impurities(V,[0,0])
bdg.set_hartree(rho)
bdg.set_fock(phi,spin_i='up',spin_f='down')
bdg.set_gorkov(chi,spin_i='up',spin_f='down')
bdg.set_hubbard_u(U,spin_i='up',spin_f='down')

bdg.solve()

hop_vector=[1,0]

# bdg.self_consistent_calculation(friction=0.7, max_iterations=1, absolute_convergence_factor=0.00001)
bdg.solve()

energy_interval=np.linspace(-4,4,51)
resolution=0.1
bdg.calculate_greens_function(energy_interval,resolution)

ldos=bdg.local_density_of_states(energy=0, atom='A')

gorkov=bdg.gorkov(atom_i=None, atom_f=None, orbital_i=None, orbital_f=None, hop_vector=[0,0], spin_i=0, spin_f=0)

ldos=np.fft.fftshift(ldos)
plt.imshow(ldos)
plt.show()
plt.close()
