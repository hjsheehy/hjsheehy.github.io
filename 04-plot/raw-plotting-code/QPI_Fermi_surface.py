from lib_plt import *
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

if len(sys.argv) < 2:
	print('Please supply name of folder')
	sys.exit()

folder = sys.argv[1]
data_location=os.path.join('../Data/',folder)

data_length = np.size(glob.glob(data_location+'\*.conf'))
# print(data_length)
def Data(i):
    global dos, omegas, index, abs_f, mu
    
    config_module = import_path(glob.glob(data_location+'\*.conf')[i]) 
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})
    
    data = np.load(glob.glob(data_location+'\*.npz')[i], allow_pickle=True)
    dos = data['dos']
    
    dos = dos.sum((2,3,4,5))
    
    index = FindNearestValueOfArray(np.real(omegas), omega)

    dos = dos[:,:,index]
    
    f = np.fft.fft2(dos, axes=(0,1), norm='ortho')
    f = np.fft.fftshift(f, axes=(0,1))
    
    abs_f = np.abs(f)
    return()
    
##########################################################
######################## Main ############################
##########################################################
def Main():
    data=[]
    mu_data=[]
    for i in range(data_length):
        Data(i)
        data.append(abs_f)
        mu_data.append(mu_xy)
    data=np.array(data)
    # max/min without central bright spot and lines:
    r=centre([0,0],n_x,n_y)
    r=(* r,)
    vmax=[]
    vmin=[]
    for i in range(data_length):
        temp=np.copy(data)
        temp[:,r[0],r[1]]=temp[:,0,0]
        vmax.append(np.amax(temp))
        vmin.append(np.amin(temp))
    vmax = np.max(vmax)
    vmin = np.max(vmin)
    print(vmax)
    print(vmin)
    vmax=vmax_set
    vmin=vmin_set
    mu_data, data = zip( *sorted( zip(mu_data, data) ) )

    # sort_index = np.argsort(mu_data)
    # mu_data=mu_data[sort_index[::-1]]
    # data=data[sort_index[::-1],:,:]

    n_rows = 10
    n_cols = 6

    data=np.reshape(data,(n_rows,n_cols,n_x,n_y))    
    mu_data=np.reshape(mu_data,(n_rows,n_cols))    
             
    xlabel='$k_x a$'
    ylabel='$k_y a$'
    
    
    cmp = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black","orange","white"])
   
    interpolation = 'none'
    extent = [-n_x/2,n_x/2,-n_y/2,n_y/2]

    fig, ax = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    for i in range(n_rows):        
        for j in range(n_cols):
            txt = f'${mu_data[i,j]:.2f}$'

            im0 = ax[i][j].imshow(
                data[i,j].T, extent=extent,
                vmin=vmin, vmax=vmax,
                interpolation=interpolation,
                cmap=cmp)    
            ax[i][j].text(0.5,0.11, 
                txt,
                {'bbox': dict(boxstyle="square", alpha=0.1, fc="white",
                               ec="none", pad=0.2)}, ha='center', va='center', c='white', transform=ax[i][j].transAxes)
            ax[i][j].set_xticklabels([])
            ax[i][j].set_yticklabels([])


    axin = inset_axes(ax[-1][-1],
                    width="10%", # width = 10% of parent_bbox width
                    height="70%", # height : 50%
                    loc=7)
    cbar2=plt.colorbar(im0, cax=axin, ticks=[vmin,vmax])
    # cbar2.ax.set_yticklabels(['Low','High'], color='white')
    cbar2.ax.set_yticklabels(['',''], color='white')

    plt.tick_params(left=False,
                    bottom=False,
                    labelleft=False,
                    labelbottom=False)
                
    title = f'QPI Fermi surfaces for various chemical potentials $\mu$'
    fig.text(0.5, 0.9,
        title, ha='center', va='center')
    fig.text(0.52, 0.09,
        r'$k_x a$', ha='center', va='center')
    fig.text(0.09, 0.5,
        r'$k_y a$', ha='center', va='center', rotation='vertical')

    fig.set_size_inches(w=latex_width, h=1.7*latex_width) 

    # plt.tight_layout()
    
    return fig

latex_width=1.2*4.7747
omega=0
vmax_set = 0.0357
vmin_set = 0

fig=Main()
plt.show()
fig.savefig('out/'+'QPI_Fermi_surface_LaNiGa2.pdf', bbox_inches = "tight")
print('Done')