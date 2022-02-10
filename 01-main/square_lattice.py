from lib import *
A=Atom([0,0],'A')
A.add_orbital('s')
lattice_vectors=[[1,0],[0,1]]
tb=Tightbinding(lattice_vectors,'Sr2RuO4')
tb.add_atom(A)
tb.n_spins=2
mu=-3.72
t=1
V=1.21
n_cells=11
tb.cut_piece(n_cells, [0,1])
tb.set_onsite(-mu,orbital='s')
tb.set_hopping(-t,hop_vector=[1,0])
tb.set_hopping(-t,hop_vector=[0,1])
tb.add_impurities(V,[0,0])

tb.solve()

energy_interval=np.linspace(-4,4,51)
resolution=0.1
tb.calculate_greens_function(energy_interval,resolution)

ldos=tb.local_density_of_states(energy=0, atom='A', orbital='s')

ldos=np.fft.fftshift(ldos)
plt.imshow(ldos)
plt.show()
plt.close()
