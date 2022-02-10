from lib import *

Ru_A=Atom([0,0],'Ru_A')
Ru_B=Atom([0.5,0],'Ru_B')
Ru_C=Atom([0.5,0.5],'Ru_C')
Ru_A.add_orbital('dxy')
Ru_A.add_orbital('dxz')
Ru_A.add_orbital('dyz')
Ru_B.add_orbital('dxy')
Ru_B.add_orbital('dxz')
Ru_C.add_orbital('dyz')
lattice_vectors=[[1,0],[0,1]]
# lattice=Crystal_lattice(lattice_vectors,'Sr2RuO4')
lattice=Tightbinding(lattice_vectors,'Sr2RuO4')
lattice.add_atom(Ru_A)
lattice.add_atom(Ru_B)
lattice.add_atom(Ru_C)
lattice.n_spins=1
lattice.cut_piece(4,0,glue_edgs=False)
# lattice.cut_piece(2,1)
mu=-2.32
# onsite=lattice._onsite_tensor(-mu,atom='Ru_B',orbital='dxz',spin=0)
so=np.eye(lattice.n_spins*Ru_B.n_orbitals)
# onsite=lattice._onsite_tensor(so,atom='Ru_B')
# print(onsite)
t=1
nn=lattice._hopping_tensor(t,k=[0.32,0.52],atom_i='Ru_A',atom_f='Ru_A',orbital_i='dxy',orbital_f='dxy',spin_i=0,spin_f=0,hop_vector=[1,-1],add_time_reversal=False)
nn=lattice._hopping_tensor(-t,k=[0.32,0.52],atom_i='Ru_A',atom_f='Ru_A',orbital_i='dxy',orbital_f='dxy',spin_i=0,spin_f=0,hop_vector=[-1,1],add_time_reversal=False)
