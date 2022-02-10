from lib_plt import *
from scipy.sparse.linalg import eigsh

if len(sys.argv) < 2:
	print('Please supply conf file')
	sys.exit()
    
data_location='../Data/Superconducting'
folder = sys.argv[1]
data_location=os.path.join(data_location,folder)
data_length = np.size(glob.glob(data_location+'\*.conf'))

def Data(i):
    global Hartree, Fock, Gorkov, U, T, mu, recursions  

    config_module = import_path(glob.glob(data_location+'\*.conf')[i]) 
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})
    
    data = np.load(glob.glob(data_location+'\*.npz')[i], allow_pickle=True)
    Hartree, Fock, Gorkov, iterations = data['Hartree'], data['Fock'], data['Gorkov'], data['iterations'] 
    return()
##########################################################
########################## Main ##########################
########################################################## 
def Main():
    Data(0)
    n_fields = 4 # phiDown, phiUp, Delta_down, Delta_up
    data = np.zeros([data_length, int(n_fields/2), int(n_fields/2), sites_away])
    r=[0,0] #since no impurity in this model
    [x0,y0]=centre(r,n_x,n_y) #centred coords
    for i in range(data_length):
        Data(i)

        data[i,0,0,:] = Gorkov[x0:x0+sites_away,y0,0,1,0,0]
        data[i,0,1,:] = -Gorkov[x0:x0+sites_away,y0,1,0,0,0]
        data[i,1,0,:] = Hartree[x0:x0+sites_away,y0,0,0]
        data[i,1,1,:] = Hartree[x0:x0+sites_away,y0,1,0]
        
    # data=Group_array(data,0,1)
    
    color=['b', 'g', 'r', 'c', 'm', 'y', 'k']
    labels=[r'$V_\sigma/t=0$', r'$V_\sigma/t=1.2$', r'$V_\downarrow/t=0$, $V_\uparrow/t=1.2$']
    marker=[ 'x', '+', '.', 'o']
    # labels=[r'No impurity', r'Impurity', r'Magnetic impurity']
    x_label = r'No. sites from $r_0$ to $r_\text{horiz}$'
    y_labels=np.array([[r'$\Delta_{{\uparrow\downarrow}}(\mathbf{{r}})$', r'$-\Delta_{{\downarrow\uparrow}}(\mathbf{{r}})$'], [r'$\phi_\uparrow(\mathbf{{r}})$', r'$\phi_\downarrow(\mathbf{{r}})$']])
    
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

    
    # ax2.annotate(r"$\Delta_0$", xy=(0, Delta0), xytext=(0.04, Delta0-0.04),
            # arrowprops=dict(arrowstyle="->"))
    # ax2.annotate(r"$T_c^\text{BCS}$", xy=(Tc, 0), xytext=(Tc+0.04, 0.06),
            # arrowprops=dict(arrowstyle="->"))
    
    # ax2.legend(title=f'$U/t$', loc="upper right", ncol = 2, fancybox=True, shadow=True)
    # ax2.set_ylabel('$\Delta/t$')
    # ax3.set_ylabel(r'$\phi_{{\uparrow,\downarrow}}/t$')
    # ax1.set_ylabel(r'Recursions')
    # ax3.set_xlabel('$k_BT/t$')

    fig.set_size_inches(w=latex_width, h=6.5) 
    return fig
latex_width=4.7747  
sites_away = 5
fig = Main()
plt.show()
# fig.savefig('out/'+folder+'.pdf', bbox_inches = "tight")