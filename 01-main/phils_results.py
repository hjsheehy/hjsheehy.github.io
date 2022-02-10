from lib import * 
#########################################################
################# Simple square bdg ###################
#########################################################
A=Atom([0,0],'A')
A.add_orbital('s')
lattice_vectors=[[1,0],[0,1]]
bdg=BogoliubovdeGennes(lattice_vectors,'Conventional BCS')
bdg.add_atom(A)
bdg.n_spins=2
mu=-1.7
V=2.21
V=0
t=1
phi=2.4
U=2.8
Delta=1.0

n_cells=11

bdg.cut_piece(n_cells, [0])
bdg.set_onsite(-mu,orbital='s')
bdg.set_hopping(-t,hop_vector=[1,0],label='t')
bdg.add_impurities(V,[0,0],label='V')
#########################################################
temperature=0.01
max_iterations=1000
friction=0.7
absolute_convergence_factor=0.00001

bdg.set_temperature(temperature)

bdg.set_hubbard_u(U*Pauli_x)
bdg.set_hartree([1.25*phi,0.75*phi])
# bdg.set_hartree(phi)
bdg.set_gorkov(Delta*1.0j*Pauli_y)

bdg.record_hartree(position_coordinates=[0,0], atom=0, orbital=0, spin='up', _print=False)
bdg.record_hartree(position_coordinates=[0,0], atom=0, orbital=0, spin='dn', _print=False)
bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i=0, atom_f=0, orbital_i=0, orbital_f=0, spin_i=0, spin_f=1, _print=False)

bdg.self_consistent_calculation(friction=friction, max_iterations=max_iterations, absolute_convergence_factor=absolute_convergence_factor)

print(bdg.free_energy[-1])

y0=np.real(bdg._hartree_iterations[0])
y1=np.real(bdg._hartree_iterations[1])
y2=np.real(bdg._gorkov_iterations[0])
y3=bdg.free_energy
plt.plot(y0,label=r'Hartree $\uparrow$')
plt.plot(y1,label=r'Hartree $\downarrow$')
plt.plot(y2,label='Gorkov')
plt.plot(y3,label='Free energy')
plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
plt.legend()
plt.show()
