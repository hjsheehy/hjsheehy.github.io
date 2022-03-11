from lib import *
V=1.28
t=1
A=Atom([0,0],'A')
A=Atom([0,0,0],'A')
A.add_orbital('s')
lattice_vectors=[[1,0],[0,1]]
lattice_vectors=[[1,0,0],[0,1,0],[0,0,1]]
# lattice_vectors=[[1]]
tb=Tightbinding(lattice_vectors,'TB')
tb.add_atom(A)
tb.n_spins=1
mu=-3.75

energy_interval=np.linspace(-2,2,401)
resolution=0.1

# n_cells=3
# tb.cut_piece(n_cells, [0], glue_edgs=True)

pts=3
x_pts=np.linspace(-np.pi, np.pi, pts)
y_pts=np.linspace(-np.pi, np.pi, pts)
z_pts=np.linspace(-np.pi, np.pi, pts)

# tb.kpts = np.array(np.meshgrid(x_pts))
# tb.kpts = np.array(np.meshgrid(x_pts,y_pts))
tb.kpts = np.array(np.meshgrid(x_pts,y_pts,z_pts))
# print(tb.kpts)
# print('#######')
# pts=pts*1j
# tb.kpts = np.array(np.mgrid[-np.pi:np.pi:pts,-np.pi:np.pi:pts,-np.pi:np.pi:pts])
# print(tb.kpts)
# exit()

# tb.set_onsite(-mu,orbital='s')
# tb.set_hopping(-t,hop_vector=[1],atom_i='A',atom_f='A', label='t')
# tb.set_hopping( lambda k : -mu-2*t*np.cos(k) )
tb.set_hopping(lambda k : -mu-2*t*np.sum(np.cos(k)))

# tb.add_impurities(V,[0],label='V')

#########################################################

tb.solve()

green = tb.calculate_greens_function(energy_interval,resolution)

##################################################

fig, axs = plt.subplots(1, 1)

########### 1D ##########

# dos = tb.density_of_states(sites='resolved', atom=0, orbital=0, spin='integrated', energy='resolved')
# fig, axs = tb.plot_band_structure_path(fig, axs, dos)

# plt.show()
# exit()

########### 2D ##########

# dos = tb.density_of_states(sites='resolved', atom=0, orbital=0, spin='integrated', energy=0)
# fig, axs = tb.plot_band_structure_2D(fig, axs, dos)

# plt.show()
# exit()

########### 3D ##########

dos = tb.density_of_states(sites='resolved', atom=0, orbital=0, spin='integrated', energy=0)
fig, axs = tb.plot_band_structure_3D(fig, axs, dos)

plt.show()
