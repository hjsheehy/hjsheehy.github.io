from lib import *
V=1.28
t=1
# A=Atom([0],'A') #1d
A=Atom([0,0],'A') #2d
B=Atom([0,0],'B') #2d
# A=Atom([0,0,0],'A') #3d
A.add_orbital('s')
B.add_orbital('s')
# lattice_vectors=[[1]] #1d
lattice_vectors=[[1,0],[0,1]] #2d
# lattice_vectors=[[1,0,0],[0,1,0],[0,0,1]] #3d
tb=Tightbinding(lattice_vectors,'TB')
tb.add_atom(A)
# tb.add_atom(B)
tb.n_spins=1
mu=2.7
s=0.

energy_interval=np.linspace(-2,2,401)
resolution=0.1

n_cells=21
tb.cut_piece(n_cells, [0,1], glue_edgs=True)

pts=11
x_pts=np.linspace(-np.pi, np.pi, pts)
y_pts=np.linspace(-np.pi, np.pi, pts)
z_pts=np.linspace(-np.pi, np.pi, pts)

# tb.kpts = np.array(np.meshgrid(x_pts)) #1d
# tb.kpts = np.array(np.meshgrid(x_pts,y_pts)) #2d
# tb.kpts = np.array(np.meshgrid(x_pts,y_pts,z_pts)) #3d

# tb.set_hopping(lambda k : -mu-2*t*np.sum(np.cos(k))-s, atom_i=0, atom_f=0)
# tb.set_hopping(lambda k : -mu-2*t*np.sum(np.cos(k))+s, atom_i=1, atom_f=1)
tb.set_onsite(-mu)
tb.set_hopping(-t,hop_vector=[1,0])
tb.set_hopping(-t,hop_vector=[0,1])

# tb.add_impurities(V,[0,0],label='V')

# tb.bulk_calculation=True
tb.bulk_calculation=False

tb.fourier_transform()
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

dos = tb.density_of_states(sites='resolved', atom='integrated', orbital='integrated', spin='integrated', energy=0)
dos = np.fft.fftshift(dos)
print(np.max(dos))
print(np.min(dos))
vmin=0
vmax=3.1
plt.imshow(dos,origin='lower',interpolation=None)#,vmin=vmin, vmax=vmax)
plt.show()
exit()
fig, axs = tb.plot_band_structure_2D(fig, axs, dos)

plt.show()

########### 3D ##########

# dos = tb.density_of_states(sites='resolved', atom=0, orbital=0, spin='integrated', energy=0)
# fig, axs = tb.plot_band_structure_3D(fig, axs, dos)

# plt.show()
