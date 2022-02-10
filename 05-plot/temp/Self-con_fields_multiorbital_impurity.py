from lib_plt import *
from scipy.sparse.linalg import eigsh

# if len(sys.argv) < 2:
	# print('Please supply conf file')
	# sys.exit()
    
# data_location='../Data/Superconducting'
# folder = sys.argv[1]
# data_location=os.path.join(data_location,folder)
# data_length = np.size(glob.glob(data_location+'\*.conf'))

# def Data(i):
    # global Hartree, Fock, Gorkov, U, T, mu, recursions  

    # config_module = import_path(glob.glob(data_location+'\*.conf')[i]) 
    # module_dict, to_import = import_all(config_module)
    # globals().update({name: module_dict[name] for name in to_import})
    
    # data = np.load(glob.glob(data_location+'\*.npz')[i], allow_pickle=True)
    # Hartree, Fock, Gorkov, iterations = data['Hartree'], data['Fock'], data['Gorkov'], data['iterations'] 
    # return()
##########################################################
########################## Temp ##########################
##########################################################
def Conf(n_spins,n_orbitals,n,mu,s,delta,t,V,impurity_locations,
    U,phi_up_plus,phi_down_plus,phi_up_minus,phi_down_minus,Delta,
    varsigma,friction,max_iterations,eps,T,
    i,j):
   
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
    impurity_spin=impurity_spin_list[i]
    impurity_orbital=V*np.eye(n_orbitals)
    impurity_tensor=np.kron(impurity_spin,impurity_orbital)
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
    Hartree[:,:,0,0]=+phi_up_plus*np.ones([n_x,n_y])
    Hartree[:,:,1,0]=+phi_down_plus*np.ones([n_x,n_y])
    # Hartree[:,:,0,1]=-phi_up_minus*np.ones([n_x,n_y])
    # Hartree[:,:,1,1]=-phi_down_minus*np.ones([n_x,n_y])

    sym = np.array([[1+varsigma, 0], [0, 1-varsigma]])
    anti_sym = 1.0j*Pauli_y
    Gorkov_spin_orbital = [np.kron(anti_sym, sym), np.kron(sym, anti_sym)][j]
    Gorkov_spin_orbital = np.kron(anti_sym,np.eye(n_orbitals))
    Gorkov_tensor=-Delta*np.kron(np.eye(n_sites), Gorkov_spin_orbital)

    ###################################################################
    ############################ Energy ###############################
    ###################################################################
    Include_DOS=True
    epsilon=0.1
    omega_pts=401
    omega_lim=10*t 
    omegas=np.linspace(-omega_lim,omega_lim,omega_pts,endpoint=True,dtype=complex)+1.0j*epsilon
    
    ############### Warning: ##################
    globals().update(locals())
    
xx=0
yy=0
zz=0
s=0
t=1
mm=0

def Update(Hartree, Fock, Gorkov):  
    global w, v
    h0 = H0(onsite_tensor, nn_tensor, nn_tensor, impurity_tensor, impurity_locations, n_x, n_y)
    ham = H_BdG(h0, Hartree, Fock, Gorkov, index, n_x,n_y,n_sites,n_spins,n_orbitals,dof)
    # del(h0)
    
    
    w,v = la.eigh(ham, overwrite_a=True)
    
    # print(np.sum(np.split(ham,2)[0]))

    dm=DM(v, n_x,n_y,index,n_sites,n_spins,n_orbitals)
    dm=Normal_Fermi(dm, dof, w, T)
    adm=ADM(v, n_x,n_y,index,n_sites,n_spins,n_orbitals)
    adm=Anomalous_Fermi(adm, dof, w, T)

    Hartree_update, Fock_update, Gorkov_update = HartreeFockGorkov(Hubbard_U, dm, adm)


    Hartree_update = (1-friction)*Hartree_update+friction*Hartree
    Fock_update = (1-friction)*Fock_update+friction*Fock
    Gorkov_update = (1-friction)*Gorkov_update+friction*Gorkov
    
    return Hartree_update, Fock_update, Gorkov_update

def Main_(Hartree, Fock, Gorkov):
    global phiUp_list, phiDown_list, GorkovUp_list, GorkovDown_list
    phiUp_list=[Hartree[xx,yy,s,mm]]
    phiDown_list=[Hartree[xx,yy,t,mm]]
    GorkovUp_list=[Gorkov[xx,yy,s,t,mm,mm]]
    GorkovDown_list=[-Gorkov[xx,yy,t,s,mm,mm]]
    for i in range(max_iterations):
        gc.collect() # Garbage collector prevents MemError
        Hartree_update, Fock_update, Gorkov_update = Update(Hartree, Fock, Gorkov) #remove 0 later
        if (np.allclose(Hartree,Hartree_update,atol=eps) & np.allclose(Fock,Fock_update,atol=eps) & np.allclose(Gorkov,Gorkov_update,atol=eps)):
            return Hartree_update, Fock_update, Gorkov_update, i+1


        DeltaUp = Gorkov_update[xx,yy,s,t,mm,mm]
        DeltaDown = Gorkov_update[xx,yy,t,s,mm,mm]
        GorkovUp_list.append(DeltaUp)
        GorkovDown_list.append(-DeltaDown)
        phiUp_list.append(Hartree_update[xx,yy,s,mm])
        phiDown_list.append(Hartree_update[xx,yy,t,mm])


        Hartree, Fock, Gorkov = Hartree_update, Fock_update, Gorkov_update
    Hartree_update, Fock_update, Gorkov_update = float('nan')*Hartree, float('nan')*Fock, float('nan')*Gorkov
    return Hartree_update, Fock_update, Gorkov_update, max_iterations

def Data(i):
    global Hartree, Fock, Gorkov, U, T, mu, iterations, index, Hubbard_U
    Conf(n_spins,n_orbitals,n,mu,s,delta,t,V,impurity_locations,
        U,phi_up_plus,phi_down_plus,phi_up_minus,phi_down_minus,Delta,
        varsigma,friction,max_iterations,eps,T,
        i,j)
    index=Index(dof,n_x,n_y,n_sites,n_spins,n_orbitals)
    Hubbard_U = indexed(Hubbard_tensor,n_x,n_y,index,n_sites,n_spins,n_orbitals)
    if 'Gorkov_tensor' in globals():
        Gorkov = indexed(Gorkov_tensor, n_x,n_y,index,n_sites,n_spins,n_orbitals)

    Hartree, Fock, Gorkov, iterations = Main_(Hartree, Fock, Gorkov)
    Gorkov=np.real_if_close(Gorkov)
    return()
##########################################################
########################## Main ##########################
########################################################## 
def Main(j):
    Data(0)
    n_fields = 4 # phiDown, phiUp, Delta_down, Delta_up
    data = np.zeros([data_length, int(n_fields/2), int(n_fields/2), sites_away])
    r=[0,0] #since no impurity in this model
    [x0,y0]=centre(r,n_x,n_y) #centred coords
    for i in range(data_length):
        Data(i)
        # print(Gorkov[x0+sites_away,y0,0,1,0,0])
        # print(Fock[x0+sites_away,y0,0,1,0,0])
        if j==0:
            data[i,0,0,:] = Gorkov[x0:x0+sites_away,y0,0,1,0,0]
            data[i,0,1,:] = -Gorkov[x0:x0+sites_away,y0,1,0,0,0]
            # data[i,1,0,:] = Gorkov[x0:x0+sites_away,y0,0,1,1,1]
            # data[i,1,1,:] = -Gorkov[x0:x0+sites_away,y0,1,0,1,1]
        else:
            data[i,0,0,:] = Gorkov[x0:x0+sites_away,y0,0,0,0,1]
            data[i,0,1,:] = -Gorkov[x0:x0+sites_away,y0,0,0,1,0]
            data[i,1,0,:] = Gorkov[x0:x0+sites_away,y0,1,1,0,1]
            data[i,1,1,:] = -Gorkov[x0:x0+sites_away,y0,1,1,1,0]

        # data[i,0,0,:] = Hartree[x0:x0+sites_away,y0,0,0]
        # data[i,0,1,:] = Hartree[x0:x0+sites_away,y0,1,0]
        # data[i,1,0,:] = Hartree[x0:x0+sites_away,y0,0,1]
        # data[i,1,1,:] = Hartree[x0:x0+sites_away,y0,1,1]
    
    color=['b', 'g', 'r', 'c', 'm', 'y', 'k']
    labels=[r'$V_\sigma/t=0$', r'$V_\sigma/t=1.2$', r'$V_\downarrow/t=0$, $V_\uparrow/t=1.2$']
    marker=[ 'x', '+', '.', 'o']
    # labels=[r'No impurity', r'Impurity', r'Magnetic impurity']
    x_label = r'No. sites from $r_0$ to $r_\text{horiz}$'
    if j==0:
        y_labels=np.array([[r'$\Delta_{{\uparrow\downarrow}}^{{++}}(\mathbf{{r}})$', r'$-\Delta_{{\downarrow\uparrow}}^{{++}}(\mathbf{{r}})$'], [r'$\Delta_{{\uparrow\downarrow}}^{{--}}(\mathbf{{r}})$', r'$-\Delta_{{\downarrow\uparrow}}^{{--}}(\mathbf{{r}})$']])
    else:
        y_labels=np.array([[r'$\Delta_{{\uparrow\uparrow}}^{{+-}}(\mathbf{{r}})$', r'$-\Delta_{{\uparrow\uparrow}}^{{-+}}(\mathbf{{r}})$'], [r'$\Delta_{{\downarrow\downarrow}}^{{+-}}(\mathbf{{r}})$', r'$-\Delta_{{\downarrow\downarrow}}^{{-+}}(\mathbf{{r}})$']])
    # y_labels=np.array([[r'$\phi_{{\uparrow}}^{{+}}(\mathbf{{r}})$', r'$\phi_{{\downarrow}}^{{+}}(\mathbf{{r}})$'], [r'$\phi_{{\uparrow}}^{{-}}(\mathbf{{r}})$', r'$\phi_{{\downarrow}}^{{-}}(\mathbf{{r}})$']])
    
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
    for i in range(data_length):
        for jx in range(int(n_fields/2)):
            for jy in range(int(n_fields/2)):
                j=jx+jy
                axs[jx, jy].plot(data[i,jx,jy,:], color=color[i], marker=marker[i], label=(labels[i]))

                axs[jx,jy].set(ylabel=y_labels[jx,jy])


    handles, labels = axs[0,0].get_legend_handles_labels()
    ncol = int(np.ceil(len(handles)/2))

    legend = fig.legend(handles, labels, loc="lower right", ncol = 3,
    fancybox=True, shadow=True,
    bbox_to_anchor=(0.44,0.9,0.5,0.5))
    fig.text(0.35, 0.05, x_label, va='center')

    fig.set_size_inches(w=latex_width, h=6.5) 
    return fig

##########################################################
########################## Conf ##########################
########################################################## 
n_spins=2
n_orbitals=1
n=11
mu=3.67
s=0#.05
delta=0#.075
t=1
V=0 #1.21
impurity_locations=[[0,0]]
U=2.8
phi_up_plus=0.9
phi_down_plus=0.9
phi_up_minus=0.7
phi_down_minus=0.7
Delta=0.5
varsigma=0#0.35
friction=0.
max_iterations=300
eps=0.001
T=0
j=0 # 0: singlet | 1:multiorbital triplet

impurity_spin_list = [np.eye(n_spins), np.eye(n_spins), np.array([[1,0],[0,0]])]
impurity_coupling_list = [0, V, V]
impurity_name_list = ['no_impurity', 'impurity', 'mag_impurity']

state_name_list = ['multiorbital_singlet', 'non-unitary_triplet']

data_length=1#len(impurity_name_list)
#############################################################
latex_width=4.7747  
sites_away = 5
fig = Main(j)
# plt.show() 
# fig.savefig('out/'+folder+'.pdf', bbox_inches = "tight")

print(GorkovUp_list[-1])
print(phiUp_list[-1])
print(len(GorkovUp_list))
