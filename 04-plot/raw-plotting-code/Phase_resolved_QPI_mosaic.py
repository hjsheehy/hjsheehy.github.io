from lib_plt import *
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

if len(sys.argv) < 2:
	print('Please supply name of folder')
	sys.exit()

folder = sys.argv[1]
data_location=os.path.join('data',folder)

# labels=np.array([
# ['''Normal_state_43_43_0_0.5_0.75_1.00_0_0_0_0_1.21''',r'Normal'],
# ['''Singlet_43_43_0_0.5_0.75_1.00_1_1_0_0_1.21''',r'Singlet'],
# ['''Unitary_triplet_43_43_0_0.5_0.75_1.00_1_1_0_0_1.21''',r'$\boldsymbol{\eta}=\mathbf{\hat{x}}$'],
# ['''Sudeep_INT_state_0.025''',r'Non-unitary triplet']])

# labels=np.array([
# ['''Normal_state_43_43_0_0.05_0.075_1.00_0_0_0_0_1.21''',r'Normal'],
# ['''Singlet_43_43_0_0.05_0.075_1.00_1_1_0_0_1.21''',r'Singlet'],
# ['''Two_orbital_Singlet_43_43_0_0.05_0.075_1.00_1_0.35_1.21''',r'Singlet $\Delta_+\not=\Delta_-$'],
# ['''Sudeep_INT_state_0.025''',r'Non-unitary triplet']])

labels=np.array([
['''Two_orbital_Singlet_43_43_0_0.05_0.075_1.00_1_0.35_1.21''',r'Multiorbital singlet, $\rho_+$',1],
['''Non-unitary_triplet_43_43_0_0.05_0.075_1.00_1_-0.35_-1j_0_1.21''',r'Non-unitary triplet, $\rho_+$',1],
['''Two_orbital_Singlet_43_43_0_0.05_0.075_1.00_1_0.35_1.21''',r'Multiorbital singlet, $\rho_-$',-1],
['''Non-unitary_triplet_43_43_0_0.05_0.075_1.00_1_-0.35_-1j_0_1.21''',r'Non-unitary triplet, $\rho_-$',-1]])

data_length=len(labels)


def Data(i):
    global data, data1, omega, dos_map, n_x, n_y
    confname=os.path.join(data_location,labels[i,0])
    confname_DOS=os.path.join(data_location,'DOS_'+labels[i,0])
    
    config_module = import_path(confname+'.conf') 
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})

    index = FindNearestValueOfArray(np.real(omegas), omega_input)
    omega=omegas[index]

    data = np.load(confname_DOS+'.npz', allow_pickle=True)
    dos = data['dos']
                
    sign = float(labels[i,2])
    ft=PhaseResolved(sign, dos, omegas, omega)
        
    dos_map=ft.sum((2,3))
    return    
##########################################################
######################## Main ############################
##########################################################
def Main():
    global data, data1, omega, dos_map
    Data(0)
    # max/min without central bright spot and lines:
    r=centre([0,0],n_x,n_y)
    r=(* r,)
    vmax=[]
    vmin=[]
    for i in range(data_length):
        Data(i)
        temp=np.copy(dos_map)
        temp[r]=temp[0,0]
        # temp[r[0],:]=temp[0,0]
        # temp[:,r[1]]=temp[0,0]
        vmax.append(np.amax(temp))
        vmin.append(np.amin(temp))
    vmax=np.amax(vmax)
    vmin=np.amin(vmin)
    print(vmax)
    print(vmin)
    # vmax=vmax_set
    # vmin=vmin_set
                
    xlabel='$k_x a$'
    ylabel='$k_y a$'
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    # newcmp = ListedColormap(cm.afmhot(np.linspace(0, 0.7, 256)))
    newcmp = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black","blue","green","yellow","red"])
    
    eV=np.real(omega)
    epsilon=np.imag(omega)
    
    text_DOS = (r'$\rho_\pm(\omega)$'
    '\n'
    f'$\omega={eV:.2f}$'
    '\n'
    f'$\epsilon={epsilon}$'
    )    
    interpolation = 'none'
    #
    extent = [-n_x/2,n_x/2,-n_y/2,n_y/2]
    
    i=0
    Data(i)
    im=ax1.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)    
    ax1.text(0.6,0.04, 
        labels[i,1],
             {'bbox': dict(boxstyle="square", alpha=0.6, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax1.transAxes)
    

    i=1
    Data(i)    
    ax2 = fig.add_subplot(222, sharey=ax1)
    ax2.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)
    ax2.text(0.6,0.04, 
        labels[i,1],
             {'bbox': dict(boxstyle="square", alpha=0.6, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax2.transAxes)


    i=2
    Data(i)
    ax3 = fig.add_subplot(223, sharex=ax1, sharey=ax1) 
    ax3.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)
    axins = inset_axes(ax3,
                    width="5%", # width = 10% of parent_bbox width
                    height="80%", # height : 50%
                    loc=6)
    cbar=plt.colorbar(im, cax=axins, ticks=[vmin,vmax])
    # vmin,vmax=round(vmin,2),round(vmax,2)
    # cbar.ax.set_yticklabels([vmin,vmax], color='white')
    cbar.ax.set_yticklabels(['Low','High'], color='white')
    ax3.text(0.6,0.04, 
        labels[i,1],
             {'bbox': dict(boxstyle="square", alpha=0.6, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax3.transAxes)

    i=3
    Data(i)
    ax4 = fig.add_subplot(224, sharex=ax2, sharey=ax1)
    im=ax4.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)
    ax4.text(0.6,0.04, 
        labels[i,1],
             {'bbox': dict(boxstyle="square", alpha=0.6, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax4.transAxes)

    ax3.set(xlabel=xlabel)
    ax4.set(xlabel=xlabel)
    ax1.set(ylabel=ylabel)
    ax3.set(ylabel=ylabel)
    
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    
    ax3.set_xticks([-n_x/2,0,n_x/2])
    ax3.set_xticklabels(['$-\pi$',0,'$\pi$'])
    ax3.set_yticks([-n_y/2,0,n_y/2])
    ax3.set_yticklabels(['$-\pi$',0,'$\pi$'])
    ax4.set_xticks([-n_x/2,0,n_x/2])
    ax4.set_xticklabels(['$-\pi$',0,'$\pi$'])
    ax1.set_yticks([-n_y/2,0,n_y/2])
    ax1.set_yticklabels(['$-\pi$',0,'$\pi$'])
    
    ax1.text(0.16,0.87, 
        text_DOS,
             {'bbox': dict(boxstyle="square", alpha=0.6, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax1.transAxes)

    fig.set_size_inches(w=latex_width, h=latex_width) 

    plt.tight_layout()
    
    print('Done')
    return fig

latex_width=4.7747
omega_input=0.6

# vmax_set=0.05
# vmin_set=-vmax_set
fig=Main()
# plt.show()
fig.savefig('out/'+f'Phase_resolved_QPI_mosaic_{np.real(omega):.2f}_'+folder+'.pdf', bbox_inches = "tight")