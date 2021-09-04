from lib_plt import *

if len(sys.argv) < 2:
	print('Please supply conf file')
	sys.exit()
    
# out_folder='../out/'
out_folder=''
    
conf = sys.argv[1]
confname = conf.split('.conf')[0]
config_module = import_path(conf)
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

def Update(Hartree, Fock, Gorkov):  
    h0 = H0(onsite_spin, onsite_orbital, nn_spin, nn_orbital, im_spin, im_orbital,
            impurity_locations, n_x, n_y)
    ham=H_BdG(h0, Hartree, Fock, Gorkov, n_x,n_y,n_sites,n_spins,n_orbitals,dof)
    w,v = la.eigh(ham)

    dm=DM(v,n_x,n_y,dof,n_sites,n_spins,n_orbitals)
    dm=Normal_Fermi(dm, dof, w, T)
    adm=ADM(v,n_x,n_y,dof,n_sites,n_spins,n_orbitals)
    adm=Anomalous_Fermi(adm, dof, w, T)

    Hartree_update, Fock_update, Gorkov_update = HartreeFockGorkov(Hubbard_U, dm, adm)

    Hartree_update = (1-f)*Hartree_update+f*Hartree
    Fock_update = (1-f)*Fock_update+f*Fock
    Gorkov_update = (1-f)*Gorkov_update+f*Gorkov
    return Hartree_update, Fock_update, Gorkov_update
##########################################################
########################## Main ##########################
########################################################## 
def Main(Hartree, Fock, Gorkov):
    Hartree_list=[]
    Fock_list=[]
    Gorkov_list=[]
    for i in range(max_data_length):
        Hartree_list.append(Hartree) 
        Fock_list.append(Fock)
        Gorkov_list.append(Gorkov)
        Hartree, Fock, Gorkov = Update(Hartree, Fock, Gorkov) #remove 0 later
                
        if (np.allclose(Hartree,Hartree_list[i],atol=eps) & np.allclose(Fock,Fock_list[i],atol=eps) & np.allclose(Gorkov,Gorkov_list[i],atol=eps)):
            return Hartree_list, Fock_list, Gorkov_list
    i=-1
    # print(Hartree[0,0,0,0])
    # x=Hartree_list[i]
    # print(np.allclose(Hartree,Hartree_list[i],eps))
    # print(x[0,0,0,0])
    print(Gorkov[:,:,:,:,0,0])
    x=Gorkov_list[i]
    print(x[0,0,:,:,0,0])
    print(np.allclose(Gorkov,x,eps))
    # print(Fock[0,0,:,:,0,0])
    # x=Fock_list[i]
    # print(x[0,0,0,0,0,0])
    # print(np.allclose(Fock,Fock_list[i],eps))
    Hartree, Fock, Gorkov = float('nan')*Hartree, float('nan')*Fock, float('nan')*Gorkov
    return Hartree_list, Fock_list, Gorkov_list


n_U=6
n_T=40
U_list = np.linspace(1.2, 2.8, n_U)
T_list = np.linspace(0, 0.5, n_T)

Hartree_data=np.zeros([n_T, n_U,n_x,n_y,n_spins,n_orbitals])
Fock_data=np.zeros([n_T, n_U,n_x,n_y,n_spins,n_spins,n_orbitals,n_orbitals])
Gorkov_data=np.zeros([n_T, n_U,n_x,n_y,n_spins,n_spins,n_orbitals,n_orbitals])

while True:
    i=0
    for T in T_list:
        j=0
        for U in U_list:
            Hubbard_tensor = -U*np.kron(np.eye(n_sites), np.kron( U_spin, U_orbital))
            Hubbard_U = indexed(Hubbard_tensor,n_x,n_y,dof,n_sites,n_spins,n_orbitals)
            
            Hartree_list, Fock_list, Gorkov_list = Main(Hartree, Fock, Gorkov)
                        
            Hartree_data[i,j,:,:,:,:]=Hartree_list[-1]
            Fock_data[i,j,:,:,:,:,:,:]=Fock_list[-1]
            Gorkov_data[i,j,:,:,:,:,:,:]=Gorkov_list[-1]
            
            j+=1
            print(len(Hartree_list))
        i+=1
        print(f"{100*i/len(T_list)}% done")
    with open(out_folder+confname+'.npz', 'wb') as f:
        np.savez(f,Hartree_data=Hartree_data,Fock_data=Fock_data,Gorkov_data=Gorkov_data,U_list=U_list,T_list=T_list)
    # test:
    ret  = zp.ZipFile(out_folder+confname+'.npz').testzip()
    if ret is not None:
        print("Broken zip: restarting calculation")
        sys.exit(1)
    else:
        break

# ax = plt.figure().gca()
# plt.plot(x,y1,color='orange',marker='+')
# plt.plot(x,y2,color='green',marker='x')
# plt.plot(x,y3,color='blue',marker='.')
# plt.plot(x,y4,color='red',marker='o')
# plt.xlim([0,len(x)])
# plt.ylim([0,3.5])
# ax.xaxis.set_major_locator(MaxNLocator(integer=True))
# plt.show()