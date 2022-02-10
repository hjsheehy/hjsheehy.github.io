import sys
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import os


"""""""""
Three band model of the Surface layer of Sr2RuO4 with spin orbit coupling
12x12 matrix
This model is taken from Lee et. al. PRB 81 1840403 2010 https://journals.aps.org/prb/pdf/10.1103/PhysRevB.81.184403
It was developed for bulk Sr3Ru2O7, however it is essentially the same as the model developed by Scaffidi et. al. (PRB 89, 220510(R) 2014 https://arxiv.org/pdf/1401.0016.pdf) except with the inclusion of the octahedral rotation.
I have modified the hopping parameters to agree with ARPES data. specifically the Fermi surface by Tamai et. al. and the bulk bands from Emil's thesis. (don't forget to put references).
"""""""""
muB = 5.788381e-5
class Sr2RuO4():
    def __init__(self,KXY_GRID,MAG_FIELD=0,NEMATIC=0,CDW=0):
        self.kxy_grid = KXY_GRID
        self.increments = np.pi/(self.kxy_grid*0.5)
        self.Nrepeats=4

        self.g = 3 #g-factor of 214
        self.Field = MAG_FIELD
        self.B = self.Field*muB*self.g/2.0 #Magnetic field
        self.Nem = NEMATIC #d-wave nematic order parameter strength

        self.CDW = CDW

        self.HamSize = 3
        self.MatSize = self.HamSize*4

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
        
        # self.mu=0
        # self.mu_c=0
        # self.eta=0
        # self.t2=0
        # self.t3=0
        # self.t4=0
        # self.t5=0
        # self.t_inter=0

    def Load_Ru_Hamiltonian(self,kx,ky):

        H_up = np.zeros((int(self.HamSize),int(self.HamSize)),dtype=np.complex)
        H_down = np.zeros((int(self.HamSize),int(self.HamSize)),dtype=np.complex)
        H = np.zeros((self.HamSize*2,self.HamSize*2),dtype=np.complex)

        E_dxz = -2*self.t1*np.cos(kx) - 2*self.t2*np.cos(ky) - self.mu
        E_dyz = -2*self.t2*np.cos(kx) - 2*self.t1*np.cos(ky) - self.mu
        E_dxy = -2*self.t3*(np.cos(kx) + np.cos(ky)) -4*self.t4*np.cos(kx)*np.cos(ky) - self.mu_c
        E_dxy += -2*self.t5*(np.cos(2*kx) +np.cos(2*ky)) #This is an additional term for the surface
        g = -4*self.t_inter*np.sin(kx)*np.sin(ky)

        #Nematicity
        E_dxz += self.Nem * (np.cos(kx) - np.cos(ky))
        E_dyz += self.Nem * (np.cos(kx) - np.cos(ky))
        E_dxy += self.Nem * (np.cos(kx) - np.cos(ky))  #d-wave nematicity

        s = 1 #spin up
        H_up[0,0] = E_dxz
        H_up[1,1] = E_dyz
        H_up[2,2] = E_dxy
        H_up[0,1] = g - 1j*self.eta*s
        H_up[0,2] = 1j*self.eta
        H_up[1,0] = g + 1j*self.eta*s
        H_up[1,2] = -self.eta*s
        H_up[2,0] = -1j*self.eta
        H_up[2,1] = -self.eta*s

        s = -1 #spin down
        H_down[0,0] = E_dxz
        H_down[1,1] = E_dyz
        H_down[2,2] = E_dxy
        H_down[0,1] = g - 1j*self.eta*s
        H_down[0,2] = 1j*self.eta
        H_down[1,0] = g + 1j*self.eta*s
        H_down[1,2] = -self.eta*s
        H_down[2,0] = -1j*self.eta
        H_down[2,1] = -self.eta*s

        #Stick the two spin block hamiltonians together.
        for i in range(int(self.HamSize)):
            for j in range(int(self.HamSize)):
                H[i,j] = H_up[i,j]
                H[i+int(self.HamSize),j+int(self.HamSize)] = H_down[i,j]

        return H


    def Load_Hamiltonian(self,kx,ky):
        k1 = kx
        k2 = ky
        kx = (k1-k2)/2
        ky = (k1+k2)/2

        """""""""
        This model is taken from Lee et. al. PRB 81 1840403 2010 https://journals.aps.org/prb/pdf/10.1103/PhysRevB.81.184403
        """""""""

        H_Surf = np.zeros((self.MatSize,self.MatSize),dtype=np.complex)
        #The two diagonal block Hamiltonians
        H_Ru1 = self.Load_Ru_Hamiltonian(kx,ky)
        H_Ru2 = self.Load_Ru_Hamiltonian(kx+np.pi,ky+np.pi)

        #Off diagonal block hamiltonian
        G = np.zeros((self.HamSize,self.HamSize),dtype=np.complex)

        G[0][0] = self.CDW
        G[1][1] = self.CDW
        G[2][2] = self.CDW

        #Add the Hamiltonians together into a full
        M = self.HamSize
        for i in range(M):
            for j in range(M):
                H_Surf[i][j] = H_Ru1[i][j] #Spin up
                H_Surf[i+M][j+M] = H_Ru2[i][j] #Spin up
                H_Surf[i+2*M][j+2*M] = H_Ru1[i+M][j+M] #Spin down
                H_Surf[i+3*M][j+3*M] = H_Ru2[i+M][j+M] #Spin down

                #Add in hybridisation (self.CDW)
                H_Surf[i+M][j] = -G[j][i] #Transpose
                H_Surf[i][j+M] = -G[i][j]
                H_Surf[i+3*M][j+2*M] = G[j][i] #Transpose
                H_Surf[i+2*M][j+3*M] = G[i][j]

        #Mag field
        for k in range(M):
            H_Surf[k][k]+= self.B
            H_Surf[M+k][M+k]+= self.B
            H_Surf[2*M+k][2*M+k]-= self.B
            H_Surf[3*M+k][3*M+k]-= self.B

        return H_Surf


    def Diagonalise(self, Hamiltonian):
        # Diagonalises Hamiltonian and returns eigenvalues, eigenvectors, the maximum orbital character, and the complete orbital character list.
        # This command needs to be updated for each individual system.
        Diag = LA.eigh(Hamiltonian)
        eigenvalues = np.around((Diag[0]), decimals=6)
        eigenvectors = Diag[1]
        store = []
        unsorted = []
        for x in range(len(eigenvalues)):
            Orb_dxz = 0
            Orb_dyz = 0
            Orb_dxy = 0
            for i in range(self.Nrepeats): #add together spin up and spin down weights on 2 Ru atoms.
                Orb_dxz += (round(abs(eigenvectors[3*i+0][x]) ** 2, 6))
                Orb_dyz += (round(abs(eigenvectors[3*i+1][x]) ** 2, 6))
                Orb_dxy += (round(abs(eigenvectors[3*i+2][x]) ** 2, 6))

            unsorted.append([Orb_dxz, Orb_dyz, Orb_dxy])
            list = [Orb_dxz, Orb_dyz, Orb_dxy]
            unsorted.append([Orb_dxz, Orb_dyz, Orb_dxy])
            list = [Orb_dxz, Orb_dyz, Orb_dxy]
            list.sort()

            if list[-1] == Orb_dxy:
                store.append("#1C397C")  # dark blue
            elif list[-1] == Orb_dxz:
                store.append("#CA0C21")  # red
            elif list[-1] == Orb_dyz:
                store.append("#058A39")  # green

        return eigenvalues,eigenvectors,store,unsorted

    def BandStructure(self,Ymin=-1,Ymax=1,mu=0):
        #Plots G-X-M-G band structure of Sr2RuO4
        kpts = [self.kxy_grid,self.kxy_grid,self.kxy_grid]

        count = 0
        EigenStore = []
        OrbStore = []
        for x in np.arange(0, np.pi, self.increments/2.0):
            count += 1
            H = self.Load_Hamiltonian(x,x)
            Diag = self.Diagonalise(H)

            Eigenvalues = Diag[0]
            OrbChar = Diag[2]
            EigenStore.append(Eigenvalues)
            OrbStore.append(OrbChar)

        for x in np.arange(0 + np.pi / self.kxy_grid, np.pi, self.increments/2.0):
            count += 1
            H = self.Load_Hamiltonian(np.pi-x,np.pi)
            Diag = self.Diagonalise(H)

            Eigenvalues = Diag[0]
            OrbChar = Diag[2]
            EigenStore.append(Eigenvalues)
            OrbStore.append(OrbChar)

        for x in np.arange(0 + np.pi / self.kxy_grid, np.pi, self.increments/2.0):
            count += 1
            H = self.Load_Hamiltonian(0,np.pi - x)
            Diag = self.Diagonalise(H)

            Eigenvalues = Diag[0]
            OrbChar = Diag[2]
            EigenStore.append(Eigenvalues)
            OrbStore.append(OrbChar)

        x = [x for x in range(len(EigenStore))]
        for i in range(len(Eigenvalues)):
            EigenVals_1 = [EigenStore[j][i] -mu for j in range(len(EigenStore))]
            OrbVals_1 = [OrbStore[j][i] for j in range(len(OrbStore))]
            plt.scatter(x, EigenVals_1, c=OrbVals_1, lw=0, s=10,zorder=1)

        plt.vlines(x=kpts[0]-1, ymin=Ymin, ymax=Ymax)
        plt.vlines(x=kpts[0]+kpts[1]-2, ymin=Ymin, ymax=Ymax)

        plt.xlim([0, kpts[0]+kpts[1]+kpts[2]-3])
        plt.ylim([Ymin, Ymax])

        plt.xticks([0, kpts[0]-1, kpts[0]+kpts[1]-2, kpts[0]+kpts[1]+kpts[2]-3], ["G", "M", "X", "G"])
        plt.ylabel("E (eV)")
        plt.xlabel("k")
        plt.axhline(0, 0, count, lw=2, color='black')
        fpdf = 'BandComp_' + str(Ymin) + "_" + str(Ymax) +".pdf"
        fpng = 'BandComp_' + str(Ymin) + "_" + str(Ymax) + ".png"
        plt.savefig(fpdf,format='pdf')
        plt.savefig(fpng,format='png',dpi=300)
        plt.show()

    def FermiSurface(self, omega, save=True):
            NBZ = 2
            Gamma = self.Load_Hamiltonian(0, 0)
            E_Gamma = self.Diagonalise(Gamma)[0]
            nr = NBZ * self.kxy_grid
            nc = NBZ * self.kxy_grid
            efermi = omega
            inc_kx = np.pi / ((nr) / (2 * NBZ))
            inc_ky = np.pi / ((nc) / (2 * NBZ))

            bands = [np.zeros((nr, nc), dtype=float) for i in range(len(E_Gamma))]

            nr = 0
            for kx in np.arange(-NBZ * np.pi, NBZ * np.pi, inc_kx):
                nc = 0
                for ky in np.arange(-NBZ * np.pi, NBZ * np.pi, inc_ky):
                    new_H = self.Load_Hamiltonian(kx, ky)
                    E_H = self.Diagonalise(new_H)[0]
                    for i in range(len(E_H)):
                        bands[i][nc][nr] = E_H[i]
                    nc += 1
                nr += 1

            for S in bands:
                x = np.arange(-NBZ * np.pi, NBZ * np.pi, inc_kx)
                y = np.arange(-NBZ * np.pi, NBZ * np.pi, inc_ky)
                X, Y = np.meshgrid(x, y)
                plt.contour(X, Y, S, [efermi], colors=[(95 / 256, 95 / 256, 95 / 256)], lw=3)


            plt.xlim([-np.pi, np.pi])
            plt.ylim([-np.pi, np.pi])
            filename = "FS"
            fpng = filename + ".png"
            fpdf = filename + ".pdf"
            if save == True:
                plt.savefig(fpng, format='png', dpi=300)
                plt.savefig(fpdf, format='pdf', dpi=64)
            plt.show()

def main():
    kxy_grid = 64 # Number of points for plotting
    Sr214 = Sr2RuO4(KXY_GRID=kxy_grid)
    Sr214.BandStructure()
    Sr214.FermiSurface(omega=0.0)

#main()

kxy_grid = 64 # Number of points for plotting
Sr214 = Sr2RuO4(KXY_GRID=kxy_grid)
Sr214.B=0
kx=0.12
ky=1.12
ham=Sr214.Load_Hamiltonian(kx,ky)
# print(np.diag(ham))
print(np.real(np.around(ham,2)))
