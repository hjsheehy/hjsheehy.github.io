from lib_plt import *
from scipy.sparse.linalg import eigsh

if len(sys.argv) < 2:
	print('Please supply conf file')
	sys.exit()
    
data_location='../Data/Singlet'
folder = sys.argv[1]
data_location=os.path.join(data_location,folder)
data_length = np.size(glob.glob(data_location+'\*.conf'))

def Data(i):
    global Hartree, Fock, Gorkov, recursions, executionTime, T, U

    config_module = import_path(glob.glob(data_location+'\*.conf')[i]) 
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})
    
    data = np.load(glob.glob(data_location+'\*.npz')[i], allow_pickle=True)
    Hartree, Fock, Gorkov, recursions, executionTime = data['Hartree'], data['Fock'], data['Gorkov'], data['recursions'], data['executionTime']
    return()
##########################################################
########################## Main ##########################
########################################################## 
def Main():
    data = np.zeros([7,data_length]) # U, T, phiDown, phiUp, Delta
    for i in range(data_length):
        Data(i)
        print(recursions)

        data[0,i], data[1,i] = U, n
        data[2,i], data[3,i] = Spin_density(Hartree, n_sites, n_orbitals)
        data[4,i] = BCS_Delta(Gorkov, n_sites, n_orbitals)
        data[5,i] = recursions
        data[6,i] = executionTime/3600
        # print(Gorkov[0,0,:,:,0,0])
        # print(Gorkov[2,2,:,:,0,0])
    data=Group_array(data,0,1)
    
    i=0
    x0 = data[1,i,:]
    label_x = r'$n$'

    y0 = data[4,i,:]
    label_y0 = r'$\Delta/t$'

    y11 = data[2,i,:]
    y12 = data[3,i,:]
    label_y1= r'$\phi_{{\uparrow,\downarrow}}/t$'

    y2 = data[5,i,:]
    label_y2= r'Recursions'    
    
    y3 = data[6,i,:]
    label_y3= r'Execution time (h)'    
    
    fig = plt.figure()

    color='tab:red'
    ax0 = fig.add_subplot(211)
    ax0.plot(x0, y0,color=color,marker='.',markersize=4,label=label_y0)
    ax0.set_ylabel(label_y0, color=color)
    ax0.tick_params(axis='y', labelcolor=color)

    ax1 = ax0.twinx()

    color = 'tab:blue'
    ax1.plot(x0, y11,color=color,marker='x',markersize=4,label=label_y1)  
    ax1.plot(x0, y12,color=color,marker='+',markersize=4,label=label_y1)  
    ax1.set_ylabel(label_y1, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    color='tab:orange'
    ax2 = fig.add_subplot(212, sharex= ax0)    
    ax2.plot(x0, y2,color=color,marker='.',markersize=4,label=label_y2)
    ax2.set_ylabel(label_y2,color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_xlabel(label_x)
    
    ax3 = ax2.twinx()

    color = 'tab:green'
    ax3.plot(x0, y3,color=color,marker='x',markersize=4,label=label_y3)  
    ax3.tick_params(axis='y', labelcolor=color)
    ax3.set_ylabel(label_y3, color=color)

    nmax=np.amax(data[1,:])
    listOf_Xticks = np.arange(3, nmax+1, 4)
    plt.xticks(listOf_Xticks)
    
    plt.setp(ax0.get_xticklabels(), visible=False)

    plt.tight_layout()

    fig.set_size_inches(w=latex_width, h=4.5) 
    return fig, ax0, ax1
    
latex_width=4.7747  
fig, ax0, ax1 = Main()
plt.show()
fig.savefig('out/'+folder+'.pdf', bbox_inches = "tight")