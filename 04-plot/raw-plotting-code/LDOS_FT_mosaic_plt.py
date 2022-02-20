from lib_plt import *
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

if len(sys.argv) < 2:
	print('Please supply name of folder')
	sys.exit()

folder = sys.argv[1]
data_location=os.path.join('data',folder)

data_length = np.size(glob.glob(data_location+'\*.conf'))

def Data():
    global data, data1, omega, dos_map
    config_module = import_path(glob.glob(data_location+'\*.conf')[0]) 
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})

    index = FindNearestValueOfArray(np.real(omegas), omega)
    omega=omegas[index]

    dos_list=[]
    for i in range(data_length):
        data = np.load(glob.glob(data_location+'\*.npz')[i], allow_pickle=True)
        dos = data['dos']
                
        dos = dos[:,:,:,:,index].sum((2,3))
        ft=FT(dos)
        
        dos_list.append(ft)
        print(f'{100*(i+1)/data_length:.2f}'+'%')

    dos_map=dos_list
    # data = np.load(os.path.join(data_location, 'Av-FT-LDOS.npz'), allow_pickle=True)
    # dos_map, omega = data['data'], data['omega']
        
    data_location1=os.path.join(data_location,'single_impurity')
    data_location1=glob.glob(data_location1+'\*.npz')[0]
    config_module = import_path(glob.glob(data_location+'\*.conf')[0]) 
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})

    data1 = np.load(data_location1, allow_pickle=True)
    dos = data1['dos']
        
    data1=FT(dos[:,:,:,:,index].sum((2,3)))
    return    
##########################################################
######################## Main ############################
##########################################################
def Main():
    global data, data1, omega, dos_map

    Data()

    k_Friedel=2*Fermi_vector(mu, t, np.real(omega))
    print(r'$k_{\text{Friedel}}$='+str(k_Friedel))

    data_length=np.shape(dos_map)[0]
    
    text1 = (f'Single central impurity')    
    text2 = (f'Trials$={1}$')    
    text3 = (f'Trials$={N}$')    
    text4 = (f'Trials$={data_length}$')    
            
    i=0
    data2=np.copy(dos_map[i])

    data3=np.copy(0*data2)
    for i in range(N):
        data3+=dos_map[i]
    data3=data3/(N)

    data4=np.copy(0*data2)
    for i in range(data_length):
        data4+=dos_map[i]
    data4=data4/(data_length)
    
    # max/min without central bright spot and lines:
    r=centre([0,0],n_x,n_y)
    r=(* r,)
    vmax=[]
    vmin=[]
    for temp_data in [data1,data2,data3,data4]:
        temp=temp_data
        temp[r]=temp[0,0]
        # vmax.append(np.amax(temp))
        # vmin.append(np.amin(temp))
    # vmax=np.amax(vmax)
    # vmin=np.amin(vmin)
    vmin=np.amin(data4)
    vmax=np.amax(data4)
    vmax1=np.amax(data1)
    vmin1=np.amin(data1)
        
    xlabel='$k_x a$'
    ylabel='$k_y a$'
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    # newcmp = ListedColormap(cm.afmhot(np.linspace(0, 0.7, 256)))
    newcmp = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black","orange","white"])
    newcmp2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black","red","white"])
    
    eV=np.real(omega)
    epsilon=np.imag(omega)
    
    text_DOS = ('$|$FT(LDOS)$|$'
    '\n'
    f'$\omega={eV:.0f}$, '
    #'\n'
    f'$\epsilon={epsilon}$'
    )    
    interpolation = 'none'
    #
    extent = [-n_x/2,n_x/2,-n_y/2,n_y/2]
    
    dos_map=data1
    im1=ax1.imshow(
        dos_map.T, extent=extent,
        vmin=vmin1, vmax=vmax1,
        interpolation=interpolation,
        cmap=newcmp2)    
    ax1.text(0.61,0.05, 
        text1,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax1.transAxes)
    axins1 = inset_axes(ax1,
                    width="5%", # width = 10% of parent_bbox width
                    height="55%", # height : 50%
                    loc=6)
    cbar1=fig.colorbar(im1, cax=axins1, ticks=[vmin1,vmax1], ax=ax1)
    cbar1.ax.set_yticklabels([f'{vmin1:.1e}',f'{vmax1:.1e}'], color='white')
    # cbar2.ax.set_yticklabels(['Low','High'], color='white')
    
    # ax1.annotate('', xy=(np.pi*(r[0]-0.5),np.pi*(r[1]+0.5)),# color='white',
            # xycoords='axes points', xytext=(k_Friedel*(r[0]+0.5), 0),
            # textcoords='offset points',
            # arrowprops=dict(arrowstyle="<-",
                            # linewidth = 2,
                            # color = 'cyan')
            # )
            
    dos_map=data2
    ax2 = fig.add_subplot(222, sharey=ax1)
    ax2.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)
    ax2.text(0.84,0.05, 
        text2,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax2.transAxes)

    dos_map=data3
    ax3 = fig.add_subplot(223, sharex=ax1, sharey=ax1) 
    im3=ax3.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)
    axins3 = inset_axes(ax3,
                    width="5%", # width = 10% of parent_bbox width
                    height="80%", # height : 50%
                    loc=6)
    cbar3=fig.colorbar(im3, cax=axins3, ticks=[vmin,vmax], ax=ax3)
    cbar3.ax.set_yticklabels([f'{vmin:.1e}',f'{vmax:.1e}'], color='white')
    # cbar.ax.set_yticklabels(['Low','High'], color='white')
    ax3.text(0.84,0.05, 
        text3,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax3.transAxes)

    dos_map=data4
    ax4 = fig.add_subplot(224, sharex=ax2, sharey=ax1)
    ax4.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)
    ax4.text(0.80,0.05, 
        text4,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
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
    
    ax1.text(0.25,0.91, 
        text_DOS,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax1.transAxes)

    fig.set_size_inches(w=latex_width, h=latex_width) 

    plt.tight_layout()
    
    print('Done')
    return fig

latex_width=4.7747
N=8
omega=0

fig=Main()
plt.show()
fig.savefig('out/'+'LDOS_FT_mosaic_'+folder+'.pdf', bbox_inches = "tight")