from lib_plt import *
from scipy.sparse.linalg import eigsh

if len(sys.argv) < 2:
	print('Please supply conf file')
	sys.exit()
    
conf = sys.argv[1]
config_module = import_path(conf)
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

index=Index(dof,n_x,n_y,n_sites,n_spins,n_orbitals)
Hubbard_U = indexed(Hubbard_tensor,n_x,n_y,index,n_sites,n_spins,n_orbitals)

def Update(Hartree, Fock, Gorkov):  
    global w, v
    h0 = H0(onsite_tensor, nn_tensor, impurity_tensor,
            impurity_locations, n_x, n_y)
    ham=H_BdG(h0, Hartree, Fock, Gorkov, index, n_x,n_y,n_sites,n_spins,n_orbitals,dof)
    w,v = la.eigh(ham)

    dm=DM(v,n_x,n_y,index,n_sites,n_spins,n_orbitals)
    dm=Normal_Fermi(dm, dof, w, T)
    adm=ADM(v,n_x,n_y,index,n_sites,n_spins,n_orbitals)
    adm=Anomalous_Fermi(adm, dof, w, T)

    Hartree_update, Fock_update, Gorkov_update = HartreeFockGorkov(Hubbard_U, dm, adm)

    Hartree_update = (1-friction)*Hartree_update+friction*Hartree
    Fock_update = (1-friction)*Fock_update+friction*Fock
    Gorkov_update = (1-friction)*Gorkov_update+friction*Gorkov
    return Hartree_update, Fock_update, Gorkov_update
##########################################################
########################## Main ##########################
########################################################## 
def Main(Hartree, Fock, Gorkov, phiUp, phiDown, Delta):
    x=[]
    y1=[]
    y2=[]
    y3=[]
    y4=[]
    for i in range(max_data_length):
        x.append(i)
        y1.append(phiUp) 
        y2.append(phiDown)
        y3.append(Delta)
        y4.append(Fock[0,0,0,0,0,0])
        Hartree, Fock, Gorkov = Update(Hartree, Fock, Gorkov) #remove 0 later
        phiDown, phiUp = Spin_density(Hartree, n_sites, n_orbitals)
        Delta = BCS_Delta(Gorkov, n_sites, n_orbitals)
        
        if (np.abs(Delta-y3[i])<eps) & (np.abs(phiUp-y1[i-1])<eps) & (np.abs(phiDown-y2[i-1])<eps):
            break
    # print(Hartree[0,0,0,0])
    # print(Fock[0,0,0,0,0,0])
    return x, y1, y2, y3, y4

# for U in [1.2,1.52,1.84,2.16,2.48,2.80]:
    # Hubbard_tensor = -U*np.kron(np.eye(n_sites), np.kron( U_spin, U_orbital))
    # Hubbard_U = indexed(Hubbard_tensor,n_x,n_y,dof,n_sites,n_spins,n_orbitals)
    
x, y1, y2, y3, y4 = Main(Hartree, Fock, Gorkov, phiUp, phiDown, Delta)

print(f'U={U}')
print(f'Iterations={len(y1)}')
print(f'phiUp={np.real(y1[-1])}')
print(f'phiDown={np.real(y2[-1])}')
print(f'Delta={np.real(y3[-1])}')


ax = plt.figure().gca()
plt.plot(x,y1,color='orange',marker='+')
plt.plot(x,y2,color='green',marker='x')
plt.plot(x,y3,color='blue',marker='.')
# plt.plot(x,y4,color='red',marker='o')
plt.xlim([0,len(x)])
plt.ylim([0,3.5])
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.show()