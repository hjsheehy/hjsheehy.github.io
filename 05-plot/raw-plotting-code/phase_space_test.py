from lib_plt import *
# out_folder='../out/'
out_folder=''

def Update(Hartree, Fock, Gorkov):  
    global w, v
    h0 = H0(onsite_spin, onsite_orbital, nn_spin, nn_orbital, im_spin, im_orbital,
            impurity_locations, n_x, n_y)
    ham=H_BdG(h0, Hartree, Fock, Gorkov, n_x,n_y,n_sites,n_spins,n_orbitals,dof)
    w,v = la.eigh(ham)

    dm=DM(v,n_x,n_y,dof,n_sites,n_spins,n_orbitals)
    dm=Normal_Fermi(dm, dof, w, T)
    adm=ADM(v,n_x,n_y,dof,n_sites,n_spins,n_orbitals)
    adm=Anomalous_Fermi(adm, dof, w, T)

    Hartree_update, Fock_update, Gorkov_update = HartreeFockGorkov(Hubbard_U, dm, adm)

    Hartree_update = (1-friction)*Hartree_update+friction*Hartree
    Fock_update = (1-friction)*Fock_update+friction*Fock
    Gorkov_update = (1-friction)*Gorkov_update+friction*Gorkov
    return Hartree_update, Fock_update, Gorkov_update
#
def Main(Hartree, Fock, Gorkov):
    for i in range(max_recursion):
        Hartree_update, Fock_update, Gorkov_update = Update(Hartree, Fock, Gorkov) #remove 0 later
        if (np.allclose(Hartree,Hartree_update,atol=eps) & np.allclose(Fock,Fock_update,atol=eps) & np.allclose(Gorkov,Gorkov_update,atol=eps)):
            return Hartree_update, Fock_update, Gorkov_update, i+1
        Hartree, Fock, Gorkov = Hartree_update, Fock_update, Gorkov_update
    Hartree_update, Fock_update, Gorkov_update = float('nan')*Hartree, float('nan')*Fock, float('nan')*Gorkov
    return Hartree_update, Fock_update, Gorkov_update, max_recursion
#
def Conf(
n_spins,
n_orbitals,
n,
mu,
s,
delta,
t,
V,
impurity_locations,
U,
phiUp,
phiDown,
DeltaUpDown, #New
DeltaDownUp, #New
friction,
max_recursion,
eps,
T,
Include_DOS):
    ###################################################################
    ######################### Normal state ############################
    ###################################################################
    onsite_spin=-mu*np.eye(n_spins)
    onsite_orbital=np.eye(n_orbitals)#-s*Pauli_z+delta*Pauli_x

    nn_spin=np.eye(n_spins)
    nn_orbital=-t*np.eye(n_orbitals)

    n_x=n_y=n
    n_sites = n_x * n_y
    dof = n_sites * n_spins * n_orbitals
    ###################################################################
    ########################## Impurities #############################
    ###################################################################
    im_spin=np.eye(n_spins)
    im_orbital=V*np.eye(n_orbitals)
    ###################################################################
    ######################### Interaction #############################
    ###################################################################
    U_orbital = np.eye(n_orbitals)
    U_spin = Pauli_x
    Hubbard_tensor = -U*np.kron(np.eye(n_sites), np.kron( U_spin, U_orbital))
    Hubbard_U = indexed(Hubbard_tensor,n_x,n_y,dof,n_sites,n_spins,n_orbitals)
    ###################################################################
    ######################### Mean-field ##############################
    ###################################################################
    Fock = np.zeros([n_x,n_y,n_spins,n_spins,n_orbitals,n_orbitals])

    Hartree = np.zeros([n_x,n_y,n_spins,n_orbitals])
    Hartree[:,:,0,:]=-phiUp*np.ones([n_x,n_y,n_orbitals])
    Hartree[:,:,1,:]=-phiDown*np.ones([n_x,n_y,n_orbitals])

    Gorkov= np.zeros([n_x,n_y,n_spins,n_spins,n_orbitals,n_orbitals])
    Gorkov[:,:,0,1,:,:]=+DeltaUpDown*np.ones([n_x,n_y,n_orbitals,n_orbitals])
    Gorkov[:,:,1,0,:,:]=-DeltaDownUp*np.ones([n_x,n_y,n_orbitals,n_orbitals])
    ###################################################################
    ############################ Energy ###############################
    ###################################################################
    epsilon=0.1
    omega_pts=401
    omega_lim=10*t 
    omegas=np.linspace(-omega_lim,omega_lim,omega_pts,endpoint=True,dtype=complex)+1.0j*epsilon
    ###################################################################
    ############### Warning: Don't try this at home  ##################
    ###################################################################
    globals().update(locals())
    
##########################################################
########################## Main ##########################
########################################################## 
n_pts=11
min,max=-1,1

n_spins=2
n_orbitals=1
n=5
mu=-2.67
s=0
delta=0
t=1
V=0
impurity_locations=[[0,0]]
U=2.8
phiUp=0.5
phiDown=0.4
DeltaUpDown_list=np.linspace(min,max,n_pts,endpoint=True)
DeltaDownUp_list=np.copy(DeltaUpDown_list)
friction=0.9
max_recursion=175
eps=0.001
T=0.165
Include_DOS=False

i=0
x_list=DeltaUpDown_list
y_list=DeltaDownUp_list
for y in y_list:
    i+=1
    for x in x_list:
        DeltaUpDown,DeltaDownUp=x,y
        Conf(n_spins,n_orbitals,n,mu,s,delta,t,V,impurity_locations,U,phiUp,phiDown,DeltaUpDown,DeltaDownUp,friction,max_recursion,eps,T,Include_DOS)
        Hartree_update, Fock_update, Gorkov_update, recursions = Main(Hartree, Fock, Gorkov)
        DeltaDownUp_update=np.sum(Gorkov_update[:,:,0,1,0,0])/n_sites
        DeltaUpDown_update=np.sum(Gorkov_update[:,:,1,0,0,0])/n_sites
        x_update_list=DeltaUpDown_update
        y_update_list=DeltaDownUp_update
        # print(recursions)
    print(f'{(100*i/len(y_list)):.2f}%: {DeltaDownUp:.4f} -> {DeltaDownUp_update:.4f}')
xx=np.array(x_list)
yy=np.array(y_list)
xx_update=np.array(x_update_list)
yy_update=np.array(y_update_list)
dx = xx_update - xx
dy = yy_update- yy

x,y = np.meshgrid(xx,yy)
u,v = np.meshgrid(dx,dy)

print(Gorkov_update[0,0,0,1,0,0])
print(Gorkov_update[0,0,1,0,0,0])
plt.quiver(x,y,u,v)
plt.show()