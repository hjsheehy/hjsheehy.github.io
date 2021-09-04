from lib_plt import *
from scipy.sparse.linalg import eigsh

if len(sys.argv) < 2:
	print('Please supply conf file')
	sys.exit()
    
data_location='../out/'
data_location=''
    
conf = sys.argv[1]
confname = conf.split('.conf')[0]
config_module = import_path(conf)
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

# data = np.load(os.path.join(data_location,'DOS_'+confname+'.npz'), allow_pickle=True)
data = np.load(os.path.join(data_location,confname+'.npz'), allow_pickle=True)
Hartree_data, Fock_data, Gorkov_data, U_list, T_list = data['Hartree_data'], data['Fock_data'], data['Gorkov_data'], data['U_list'], data['T_list']
##########################################################
########################## Main ##########################
########################################################## 
def Main():
    n_T=len(T_list)
    n_U=len(U_list)

    phiDown_data=np.zeros([n_T,n_U])
    phiUp_data=np.zeros([n_T,n_U])
    Delta_data=np.zeros([n_T,n_U])
    for i in range(n_T):
        for j in range(n_U):
            phiDown_data[i,j], phiUp_data[i,j] = Spin_density(Hartree_data[i,j], n_sites, n_orbitals)
            Delta_data[i,j] = BCS_Delta(Gorkov_data[i,j], n_sites, n_orbitals)

    color=['b', 'g', 'r', 'c', 'm', 'y', 'k']
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    for i in range(n_U):
        l1,=ax1.plot(T_list, Delta_data[:,i], color=color[i], label=f'${U_list[i]:.2f}$')
        l1,=ax1.plot(T_list, Delta_data[:,i], color=color[i], marker='.')
        l2,=ax2.plot(T_list, phiDown_data[:,i], color=color[i], marker='x')
        l3,=ax2.plot(T_list, phiUp_data[:,i], color=color[i], marker='+')

    ax1.legend(title=f'$U/t$', loc="upper right", ncol = 2, fancybox=True, shadow=True)
    ax1.set_ylabel('$\Delta/t$')
    ax2.set_ylabel(r'$\phi_{{\uparrow,\downarrow}}/t$')
    ax2.set_xlabel('$k_BT/t$')

    i=n_U-1
    curve=Delta_data[:,i]
    Delta0=curve[0]
    Tc=Delta0/1.76
    print(Tc)

    ax1.annotate(r"$\Delta_0$", xy=(0, Delta0), xytext=(0.04, Delta0-0.04),
            arrowprops=dict(arrowstyle="->"))
    ax1.annotate(r"$T_c^\text{BCS}$", xy=(Tc, 0), xytext=(Tc+0.04, 0.06),
            arrowprops=dict(arrowstyle="->"))

    fig.set_size_inches(w=latex_width, h=4.5) 

    plt.show()
latex_width=4.7747  
Main()