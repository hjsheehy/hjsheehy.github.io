from lib import *
A=Atom([0,0],'A')
B=Atom([0.25,1],'B')
A.add_orbital('s')
B.add_orbital('s')
lattice_vectors=[[1,0],[0,1]]
tb=Tightbinding(lattice_vectors,'Sr2RuO4')
tb.add_atom(A)
tb.add_atom(B)
tb.n_spins=2
mu=-3.72
t=1
V=10.21
n_cells=11
tb.cut_piece(n_cells, [0])
tb.set_onsite(-mu,orbital='s')
tb.set_hopping(-t,hop_vector=[0,0],atom_i='A',atom_f='B', label='v')
tb.set_hopping(-t,hop_vector=[1,0],atom_i='B',atom_f='A', label='w')
tb.add_impurities(V,[0,0],label='V')

tb.solve()

energy_interval=np.linspace(-4,4,101)
resolution=0.1
tb.calculate_greens_function(energy_interval,resolution)

fig,ax = plt.subplots(1,figsize=(LATEX_WIDTH,LATEX_WIDTH))

tb.plot_lattice(fig, ax, energy=None, atoms=None, plot_ldos=True, plot_magnetism=False, s=100)
# fig,ax = tb.plot_unit_cell(fig,ax)

plt.show()

# print(tb.impurities)
# print(tb.impurity_labels)
# print(tb.hoppings)

# tb.solve()

# energy_interval=np.linspace(-4,4,51)
# resolution=0.1
# tb.calculate_greens_function(energy_interval,resolution)

# ldos=tb.local_density_of_states(energy=0, atom='A', orbital='s')

# ldos=np.fft.fftshift(ldos)
# plt.imshow(ldos)
# plt.show()
# plt.close()
