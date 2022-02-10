from lib import *

class Sr2RuO4_TB():

    def __init__(self):

        self.tomeV = 0.15
        self.t1 = 1.0*self.tomeV #nearest neighbour dxz-dxz hopping along x-direction or dyz-dyz hopping along y-direction
        self.t2 = 0.1*self.t1    #nearest neighbour dxz-dxz hopping along y-direction or dyz-dyz hopping along x-direction
        self.t3 = 0.8*self.t1    #nearest neighbour dxy-dxy hopping along x and y direction
        self.t4 = 0.3*self.t1    #next nearest neighbour dxy-dxy hopping along x and y direction
        self.t5 = 0.095*self.t1 #Bulk value is 0*self.t1  #third nearest neighbour dxy-dxy hopping along x and y direction
        self.t_inter = 0.01*self.t1 #dxz-dyz hopping (nearest neighbour)
        self.mu = 0.75*self.t1  #Bulk value is 1.0*self.t1 Chemical potential
        self.mu_c = 0.812*self.t1 #Bulk value is 1.1*self.t1 Crystal field splitting between dxz/dyz and dxy orbital
        self.eta = 0.1*self.t1 #SOC

    def Load_Hamiltonian(self,kx,ky):
        
        # E_dxz = -2*self.t1*np.cos(kx) - 2*self.t2*np.cos(ky) - self.mu
        # E_dyz = -2*self.t2*np.cos(kx) - 2*self.t1*np.cos(ky) - self.mu
        # E_dxy = -2*self.t3*(np.cos(kx) + np.cos(ky)) -4*self.t4*np.cos(kx)*np.cos(ky) - self.mu_c
        # E_dxy += -2*self.t5*(np.cos(2*kx) +np.cos(2*ky)) #This is an additional term for the surface
        # g = -4*self.t_inter*np.sin(kx)*np.sin(ky)

        # model = tight_binding(dimensions=[1,1], n_spins=2, basis=[[1,0],[0,1]], orbitals=[[1,0],[1,0],[0,1]], pbc=[True,True])

        # model.set_onsite(E_dxz,orbital=0)
        # model.set_onsite(E_dyz,orbital=1)
        # model.set_onsite(E_dxy,orbital=2)

        # model.set_hopping([g-1j*self.eta*(-1)**s for s in range(2)], orb_i=0, orb_f=1)
        # model.set_hopping(1j*self.eta, orb_i=0, orb_f=2)
        # model.set_hopping([-self.eta*(-1)**s for s in range(2)], orb_i=1, orb_f=2)

        ################
        self.mu=0
        self.mu_c=0
        self.eta=0
        # self.t1=0
        self.t2=0
        self.t3=0
        self.t4=0
        self.t5=0
        self.t_inter=0
        ################

        k=np.array([kx,ky])

        Ru_A=Atom([0,0],'Ru_A')
        Ru_A.add_orbital('dxz')
        Ru_A.add_orbital('dyz')
        Ru_A.add_orbital('dxy')
        Ru_B=Atom([1,0],'Ru_B')
        Ru_B.add_orbital('dxz')
        Ru_B.add_orbital('dyz')
        Ru_B.add_orbital('dxy')
        lattice_vectors=[[2,0],[0,1]]
        lattice=Tightbinding(lattice_vectors,'Sr2RuO4')
        lattice.add_atom(Ru_A)
        lattice.add_atom(Ru_B)

        matr=lattice._onsite_tensor(-self.mu,orbital='dxz')
        matr+=lattice._onsite_tensor(-self.mu,orbital='dyz')
        matr+=lattice._onsite_tensor(-self.mu_c,orbital='dxy')

        k2=k+np.pi
        nn=lattice._hopping_tensor(-self.t1,k=k,atom_i='Ru_A',atom_f='Ru_B',orbital_i='dxz',orbital_f='dxz',hop_vector=[0,0])
        # nn+=lattice._hopping_tensor(-self.t1,k=k2,atom_i='Ru_B',atom_f='Ru_A',orbital_i='dxz',orbital_f='dxz',hop_vector=[1,0])
        # nn+=lattice._hopping_tensor(-self.t2,k=k,atom_i='Ru_A',atom_f='Ru_B',orbital_i='dyz',orbital_f='dyz',hop_vector=[0,0])
        # nn+=lattice._hopping_tensor(-self.t2,k=k2,atom_i='Ru_B',atom_f='Ru_A',orbital_i='dyz',orbital_f='dyz',hop_vector=[1,0])
        # nn+=lattice._hopping_tensor(-self.t1,k=k,atom_i='Ru_A',atom_f='Ru_A',orbital_i='dyz',orbital_f='dyz',hop_vector=[0,1])
        # nn+=lattice._hopping_tensor(-self.t1,k=k2,atom_i='Ru_B',atom_f='Ru_B',orbital_i='dyz',orbital_f='dyz',hop_vector=[0,1])
        nn+=lattice._hopping_tensor(-self.t2,k=k,atom_i='Ru_A',atom_f='Ru_A',orbital_i='dxz',orbital_f='dxz',hop_vector=[0,1])
        nn+=lattice._hopping_tensor(-self.t2,k=k2,atom_i='Ru_B',atom_f='Ru_B',orbital_i='dxz',orbital_f='dxz',hop_vector=[0,1])
        nn+=lattice._hopping_tensor(-self.t3,k=k,atom_i='Ru_A',atom_f='Ru_A',orbital_i='dxy',orbital_f='dxy',hop_vector=[0,1])
        nn+=lattice._hopping_tensor(-self.t3,k=k,atom_i='Ru_A',atom_f='Ru_A',orbital_i='dxy',orbital_f='dxy',hop_vector=[1,0])
        nn+=lattice._hopping_tensor(-self.t3,k=k2,atom_i='Ru_B',atom_f='Ru_B',orbital_i='dxy',orbital_f='dxy',hop_vector=[0,1])
        nn+=lattice._hopping_tensor(-self.t3,k=k2,atom_i='Ru_B',atom_f='Ru_B',orbital_i='dxy',orbital_f='dxy',hop_vector=[1,0])
        nn+=lattice._hopping_tensor(-self.t4,k=k,atom_i='Ru_A',atom_f='Ru_A',orbital_i='dxy',orbital_f='dxy',hop_vector=[1,1])
        nn+=lattice._hopping_tensor(-self.t4,k=k,atom_i='Ru_A',atom_f='Ru_A',orbital_i='dxy',orbital_f='dxy',hop_vector=[1,-1])
        nn+=lattice._hopping_tensor(-self.t4,k=k2,atom_i='Ru_B',atom_f='Ru_B',orbital_i='dxy',orbital_f='dxy',hop_vector=[1,1])
        nn+=lattice._hopping_tensor(-self.t4,k=k2,atom_i='Ru_B',atom_f='Ru_B',orbital_i='dxy',orbital_f='dxy',hop_vector=[1,-1])
        nn+=lattice._hopping_tensor(-self.t5,k=k,atom_i='Ru_A',atom_f='Ru_A',orbital_i='dxy',orbital_f='dxy',hop_vector=[2,0])
        nn+=lattice._hopping_tensor(-self.t5,k=k,atom_i='Ru_A',atom_f='Ru_A',orbital_i='dxy',orbital_f='dxy',hop_vector=[0,2])
        nn+=lattice._hopping_tensor(-self.t5,k=k2,atom_i='Ru_B',atom_f='Ru_B',orbital_i='dxy',orbital_f='dxy',hop_vector=[2,0])
        nn+=lattice._hopping_tensor(-self.t5,k=k2,atom_i='Ru_B',atom_f='Ru_B',orbital_i='dxy',orbital_f='dxy',hop_vector=[0,2])
        nn+=lattice._hopping_tensor(self.t_inter,k=k,atom_i='Ru_A',atom_f='Ru_A',orbital_i='dxz',orbital_f='dyz',hop_vector=[1,1])
        nn+=lattice._hopping_tensor(-self.t_inter,k=k,atom_i='Ru_A',atom_f='Ru_A',orbital_i='dxz',orbital_f='dyz',hop_vector=[1,-1])
        nn+=lattice._hopping_tensor(self.t_inter,k=k2,atom_i='Ru_B',atom_f='Ru_B',orbital_i='dxz',orbital_f='dyz',hop_vector=[1,1])
        nn+=lattice._hopping_tensor(-self.t_inter,k=k2,atom_i='Ru_B',atom_f='Ru_B',orbital_i='dxz',orbital_f='dyz',hop_vector=[1,-1])

        nn+=lattice._hopping_tensor([-1j*self.eta*(-1)**s for s in range(2)], atom_i='Ru_A',atom_f='Ru_A', orbital_i='dxz', orbital_f='dyz')
        nn+=lattice._hopping_tensor(1j*self.eta, atom_i='Ru_A',atom_f='Ru_A', orbital_i='dxz', orbital_f='dyz')
        nn+=lattice._hopping_tensor([-self.eta*(-1)**s for s in range(2)], atom_i='Ru_A',atom_f='Ru_A', orbital_i='dxz', orbital_f='dyz')
        nn+=lattice._hopping_tensor([-1j*self.eta*(-1)**s for s in range(2)], atom_i='Ru_B',atom_f='Ru_B', orbital_i='dxz', orbital_f='dyz')
        nn+=lattice._hopping_tensor(1j*self.eta, atom_i='Ru_B',atom_f='Ru_B', orbital_i='dxz', orbital_f='dyz')
        nn+=lattice._hopping_tensor([-self.eta*(-1)**s for s in range(2)], atom_i='Ru_B',atom_f='Ru_B', orbital_i='dxz', orbital_f='dyz')

        return matr+nn

#         return model._hamiltonian

Sr214=Sr2RuO4_TB()
kx=ky=0
ham=Sr214.Load_Hamiltonian(kx=kx,ky=ky)
# print(np.diag(ham))
print(ham)
