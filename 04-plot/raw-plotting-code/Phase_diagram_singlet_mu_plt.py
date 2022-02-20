from lib_plt import *
from scipy.sparse.linalg import eigsh

if len(sys.argv) < 2:
	print('Please supply location')
	sys.exit()
    
data_location='../Data/Superconducting'
folder = sys.argv[1]
data_location=os.path.join(data_location,folder)
data_length = np.size(glob.glob(data_location+'\*.conf'))
print(data_length)
def Data(i):
    global Hartree, Fock, Gorkov, U, T, mu, iterations,to_import
    config_module = import_path(glob.glob(data_location+'\*.conf')[i]) 
    module_dict, to_import = import_all(config_module)
    
    #del 'Hubbard_tensor' (excessive memory)
    index = to_import.index('Hubbard_tensor')
    to_import.pop(index)
    del module_dict['Hubbard_tensor']
        
    globals().update({name: module_dict[name] for name in to_import})

    data = np.load(glob.glob(data_location+'\*.npz')[i], allow_pickle=True)
    Hartree, Fock, Gorkov, iterations = data['Hartree'], data['Fock'], data['Gorkov'], data['iterations'] 
    return()
Data(0)
##########################################################
########################## Main ##########################
########################################################## 
def Main():
    data = np.zeros([6,data_length]) # U, T, phiDown, phiUp, Delta, f
    for i in range(data_length):
        Data(i)
        
        data[0,i], data[1,i] = U, mu
        # data[0,i], data[1,i] = U, T
        data[2,i], data[3,i] = Spin_density(Hartree, n_sites, n_orbitals)
        data[4,i] = np.abs(BCS_Delta(Gorkov, n_sites, n_orbitals))
        data[5,i] = iterations
 
    data=Group_array(data,0,1)
    
    color=['b', 'g', 'r', 'c', 'm', 'y', 'k']

    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    
    print(np.shape(data))
    n_U=np.shape(data)[1]
    for i in range(n_U):
        l3=ax3.plot(data[1,i,:], data[2,i,:], color=color[i], label=f'${data[0,i,0]:.2f}$')
        # l2=ax2.plot(data[1,i,:], data[4,i,:], color=color[i], label=f'${data[0,i,0]:.2f}$')
        
        l1=ax1.plot(data[1,i,:], data[5,i,:], color=color[i], marker='.')
        l2=ax2.plot(data[1,i,:], data[4,i,:], color=color[i], marker='.')
        l3=ax3.plot(data[1,i,:], data[2,i,:], color=color[i], marker='x')
        l4=ax3.plot(data[1,i,:], data[3,i,:], color=color[i], marker='+')
        
    # for i in range(n_U):
        # curve=data[4,i,:]
        # Delta0=curve[0]
        # Tc=Delta0/1.76
    # print(f'Delta={Delta0}')
    # print(f'Tc={Tc}')
    
    # ax2.annotate(r"$\Delta_0$", xy=(0, Delta0), xytext=(0.04, Delta0-0.04),
            # arrowprops=dict(arrowstyle="->"))
    # ax2.annotate(r"$T_c^\text{BCS}$", xy=(Tc, 0), xytext=(Tc+0.04, 0.06),
            # arrowprops=dict(arrowstyle="->"))
    
    ax3.legend(title=f'$U/t$', loc="upper left", ncol = 2, fancybox=True, shadow=True)
    # ax2.legend(title=f'$U/t$', loc="upper right", ncol = 2, fancybox=True, shadow=True)
    ax2.set_ylabel('$\Delta/t$')
    ax3.set_ylabel(r'$\phi_{{\uparrow,\downarrow}}/t$')
    ax1.set_ylabel(r'Iterations')
    ax3.set_xlabel('$\mu/t$')
    # ax3.set_xlabel('$k_BT/t$')

    fig.set_size_inches(w=latex_width, h=6.5) 
    return fig, (ax1, ax2)
latex_width=4.7747  
fig, (ax1, ax2) = Main()
plt.show()
fig.savefig('out/'+folder+'.pdf', bbox_inches = "tight")