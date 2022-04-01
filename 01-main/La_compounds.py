from lib import *
A=Atom([0,0],'A')
B=Atom([0.5,0],'B')
A.add_orbital('s')
B.add_orbital('s')
lattice_vectors=[[1,0],[0,1]]
bdg=BogoliubovdeGennes(lattice_vectors,'La_compounds')
bdg.add_atom(A)
bdg.add_atom(B)
bdg.n_spins=2
mu=3.72
t=1
v=1.21
w=0.78
Uv=2.28
Uw=1.88
rho=2.21
phi=0.32
chi=1.42
V=1.21
n_cells=11
bdg.cut_piece(n_cells, [0,1])
bdg.set_onsite(-mu)
bdg.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B',label='v')
bdg.set_hopping(-w,hop_vector=[1,0],atom_i='B',atom_f='A',label='w')
bdg.set_hopping(-t,hop_vector=[0,1],atom_i='A',atom_f='A',label='t')
bdg.set_hopping(-t,hop_vector=[0,1],atom_i='B',atom_f='B',label='t')
bdg.add_impurities(V,[0,0])
bdg.set_hartree(rho)
bdg.set_fock(phi,atom_i='A',atom_f='B')
bdg.set_fock(phi,atom_i='B',atom_f='A',hop_vector=[1,0])
bdg.set_gorkov(chi,atom_i='A',atom_f='B')
bdg.set_gorkov(chi,atom_i='B',atom_f='A',hop_vector=[1,0])

bdg.set_hubbard_u(-Uv,atom_i='A',atom_f='B',hop_vector=[0,0],orbital_i='s',orbital_f='s',spin_i=0,spin_f=0)
bdg.set_hubbard_u(Uw,atom_i='B',atom_f='A',hop_vector=[1,0])

bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B', orbital_i='s', orbital_f='s', spin_i=0, spin_f=0,_print=False)
bdg.record_gorkov(location_i=[0,0], location_f=[1,0], atom_i='B', atom_f='A', orbital_i='s', orbital_f='s', spin_i=0, spin_f=0,_print=False)

bdg.self_consistent_calculation(friction=0., max_iterations=1, absolute_convergence_factor=0.00001)

energy_interval=np.linspace(-4,4,51)
resolution=0.1
bdg.calculate_greens_function(energy_interval,resolution)

gorkov_v=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])
gorkov_w=bdg.gorkov(atom_i='B', atom_f='A', hop_vector=[1,0])

plt.plot(np.real(bdg._gorkov_iterations[0]))
plt.plot(np.real(bdg._gorkov_iterations[1]))
plt.show()
plt.close()

fig, ax = plt.subplots(1, 1)
bdg.plot_unit_cell(fig, ax, atoms='all', s=10)
plt.show()
plt.close()

ldos=bdg.local_density_of_states(energy=0, atom='A')
ldos=np.fft.fftshift(ldos)
plt.imshow(ldos)
plt.show()
plt.close()

ldos=bdg.local_density_of_states(energy=0, atom='B')
ldos=np.fft.fftshift(ldos)
plt.imshow(ldos)
plt.show()
plt.close()
