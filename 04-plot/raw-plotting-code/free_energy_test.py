from lib_plt import *
# out_folder='../out/'
out_folder=''

def Update(Hartree, Fock, Gorkov):  
    global w0, w, v
    h0 = H0(onsite_tensor, nn_tensor, impurity_tensor,
            impurity_locations, n_x, n_y)
    w0,v0=la.eigh(h0)
    ham= H_BdG(h0, Hartree, Fock, Gorkov, index,
            n_x,n_y,n_sites,n_spins,n_orbitals,dof)
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
#
def F_MF(w0, w, T, dof):
    Eg = np.sum(w0+w[dof:])
    if T==0:
        return Eg
    else:
        return Eg-T*np.sum(np.log(1+np.exp(-w/T)))
#
def FreeEnergy(Hartree, Fock, Gorkov, Hartree_update, Fock_update, Gorkov_update, U, F_MF,
    n_sites, n_spins, n_orbitals):
    chi    = np.sum(Gorkov_update[:,:,0,1,:,:])/(n_sites*n_orbitals)
    chiBar = np.sum(Gorkov_update[:,:,1,0,:,:])/(n_sites*n_orbitals)
    phi    = np.sum(Hartree_update[:,:,0,:])/(n_sites*n_orbitals)
    phiBar = np.sum(Hartree_update[:,:,1,:])/(n_sites*n_orbitals)
    Delta_    = np.sum(Gorkov[:,:,0,1,:,:])/(n_sites*n_orbitals)
    DeltaBar_ = -np.sum(Gorkov[:,:,1,0,:,:])/(n_sites*n_orbitals)
    N     = np.sum(Hartree[:,:,0,:])/(n_sites*n_orbitals)
    NBar  = np.sum(Hartree[:,:,1,:])/(n_sites*n_orbitals)
    return (-U*(Delta_*DeltaBar_+N*NBar)+phi*NBar+phiBar*N+chi*DeltaBar_+chiBar*Delta_ + F_MF)
#
# def Main(Hartree, Fock, Gorkov):
    # for i in range(max_recursion):
        # Hartree_update, Fock_update, Gorkov_update = Update(Hartree, Fock, Gorkov) #remove 0 later
        # if (np.allclose(Hartree,Hartree_update,atol=eps) & np.allclose(Fock,Fock_update,atol=eps) & np.allclose(Gorkov,Gorkov_update,atol=eps)):
            # return Hartree_update, Fock_update, Gorkov_update, i+1
        # Hartree, Fock, Gorkov = Hartree_update, Fock_update, Gorkov_update
    # Hartree_update, Fock_update, Gorkov_update = float('nan')*Hartree, float('nan')*Fock, float('nan')*Gorkov
    # return Hartree_update, Fock_update, Gorkov_update, max_recursion
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
Delta,
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
    onsite_tensor=np.kron(onsite_spin,onsite_orbital)    
    
    nn_spin=np.eye(n_spins)
    nn_orbital=-t*np.eye(n_orbitals)
    nn_tensor=np.kron(nn_spin,nn_orbital)

    n_x=n_y=n
    n_sites = n_x * n_y
    dof = n_sites * n_spins * n_orbitals
    ###################################################################
    ########################## Impurities #############################
    ###################################################################
    im_spin=np.eye(n_spins)
    im_orbital=V*np.eye(n_orbitals)
    impurity_tensor=np.kron(im_spin,im_orbital)
    ###################################################################
    ######################### Interaction #############################
    ###################################################################
    U_orbital = np.eye(n_orbitals)
    U_spin = Pauli_x
    Hubbard_tensor = -U*np.kron(np.eye(n_sites), np.kron( U_spin, U_orbital))
    ###################################################################
    ######################### Mean-field ##############################
    ###################################################################
    Fock = np.zeros([n_x,n_y,n_spins,n_spins,n_orbitals,n_orbitals])

    Hartree = np.zeros([n_x,n_y,n_spins,n_orbitals])
    Hartree[:,:,0,:]=-phiUp*np.ones([n_x,n_y,n_orbitals])
    Hartree[:,:,1,:]=-phiDown*np.ones([n_x,n_y,n_orbitals])

    Gorkov= np.zeros([n_x,n_y,n_spins,n_spins,n_orbitals,n_orbitals])
    Gorkov[:,:,0,1,:,:]=+Delta*np.ones([n_x,n_y,n_orbitals,n_orbitals])
    Gorkov[:,:,1,0,:,:]=-Delta*np.ones([n_x,n_y,n_orbitals,n_orbitals])
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
n_pts=21
T_pts=5
min,max=-.8,.8

n_spins=2
n_orbitals=1
n=11
mu=-0.67
s=0
delta=0
t=1
V=0
impurity_locations=[[0,0]]
U=2.8
phiUp=0.5
phiDown=0.5
Deltas=np.linspace(min,max,n_pts,endpoint=True)
friction=0.9
max_recursion=175
eps=0.001
TT=np.linspace(0,0.5,T_pts,endpoint=True)
Include_DOS=False

i=0
xx=Deltas
yy=np.zeros(np.shape(xx))
zz=np.zeros(np.shape(xx))
ww=np.zeros(np.shape(xx))
color=['red','green','blue']
# for T in TT:
T=0.05
for i in range(len(xx)):
    Delta=xx[i]
    Conf(n_spins,n_orbitals,n,mu,s,delta,t,V,impurity_locations,U,phiUp,phiDown,Delta,friction,max_recursion,eps,T,Include_DOS)
    index = Index(dof, n_x, n_y, n_sites,n_spins,n_orbitals)
    Hubbard_U = indexed(Hubbard_tensor,n_x,n_y,index,n_sites,n_spins,n_orbitals)
    Hartree_update, Fock_update, Gorkov_update = Update(Hartree, Fock, Gorkov)
    Delta_update=np.sum(Gorkov_update[:,:,0,1,0,0])/n_sites
    F_mf = F_MF(w0, w, T, dof)
    # Eg = np.sum(w0+w[:dof])
    # yy[i] = Eg
    zz[i] = F_mf
    ww[i] = FreeEnergy(Hartree, Fock, Gorkov, Hartree_update, Fock_update, Gorkov_update, U, F_mf, n_sites, n_spins, n_orbitals)
# plt.plot(xx,yy,label=r'$E_g$')
plt.plot(xx,zz,label=r'$F_{MF}$',color=color[0])
plt.plot(xx,ww-zz,label=r'$F-F_{MF}$')#f'$T={T}$')#,color=color[2])
plt.plot(xx,ww,label=r'$F$')
plt.plot(xx,ww*0,color='black')
plt.legend()
plt.show()