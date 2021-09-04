from lib_plt import *
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def Main():
    config_module = import_path(data_location+'.conf')
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})

    #################
    # index=Index(dof, n_x, n_y, n_sites,n_spins,n_orbitals)
    #####################################################################
    ############################ Main ###################################
    #####################################################################
    # h0 = H0(onsite_tensor, nn_tensor_x, nn_tensor_y, impurity_tensor,impurity_locations, n_x, n_y)
    # mf = Mean_field(SC_tensor, n_sites)
    # ham = H(h0, mf)
    # del(h0)
    # del(mf)
    # w,v = la.eigh(ham, overwrite_a=True)

    # density_matrix = DM(v, n_x,n_y,index,n_sites,n_spins,n_orbitals)
    # dos = DOS(omegas, w, density_matrix)
    ########################
   
    data = np.load(data_location+'.npz', allow_pickle=True)
    dos = data['dos']
    
    dos=np.trace(dos,axis1=2,axis2=3)
    dos=np.trace(dos,axis1=2,axis2=3)
    # m=0
    # dos=dos[:,:,m,m]
        
    index = FindNearestValueOfArray(np.real(omegas), omega)

    dos = dos[:,:,index]
    
    f = np.fft.fft2(dos, axes=(0,1), norm='ortho')
    f = np.fft.fftshift(f, axes=(0,1))
    
    abs_f = np.abs(f)

    data = abs_f
    # data = dos
    # max/min without central bright spot and lines:
    r=centre([0,0],n_x,n_y)
    r=(* r,)
    vmax=[]
    vmin=[]
    temp=np.copy(data)
    temp[r[0],r[1]]=temp[0,0]
    vmax.append(np.amax(temp))
    vmin.append(np.amin(temp))
    vmax = np.max(vmax)
    vmin = np.max(vmin)
    print(vmax)
    print(vmin)
    # vmax=vmax_set
    # vmin=vmin_set
                 
    cmp = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black","orange","white"])
   
    interpolation = 'none'
    extent = [-n_x/2,n_x/2,-n_y/2,n_y/2]

    fig, ax = plt.subplots(1, 1)
    
    im = ax.imshow(
        data.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=cmp)    

    ax.set_xticklabels([])
    ax.set_yticklabels([])

    axin = inset_axes(ax,
                    width="5%", # width = 10% of parent_bbox width
                    height="70%", # height : 50%
                    loc=6)
    cbar=plt.colorbar(im, cax=axin, ticks=[vmin,vmax])
    cbar.ax.set_yticklabels(['Low','High'], color='white')
    # cbar.ax.set_yticklabels(['',''], color='white')

    plt.tick_params(left=False,
                    bottom=False,
                    labelleft=False,
                    labelbottom=False)
                
    title = f'QPI($\omega={omega:.2f})$'
    ax.text(0.01,0.978, 
        title,
        {'bbox': dict(boxstyle="square", alpha=0.5, fc="black",
                       ec="none", pad=0.2)}, ha='left', va='center', c='white', transform=ax.transAxes)
    fig.text(0.5, 0.02,
        xlabel, ha='center', va='bottom')
    fig.text(-0.01, 0.5,
        ylabel, ha='left', va='center', rotation='vertical')

    fig.set_size_inches(w=latex_width, h=1.1*latex_width) 

    plt.tight_layout()
    
    return fig

xlabel=r'$k_x a$'
ylabel=r'$k_z c$'
latex_width=4.7747
omega=0
vmax_set = 0.008
vmin_set = 0
# data_location='../Data/LaNiGa2/LaNiGa2'
data_location='../Data/LaNiGa2/LaNiGa2_N'

fig=Main()
plt.show()
fig.savefig('out/'+'QPI_LaNiGa2.pdf', bbox_inches = "tight")
print('Done')