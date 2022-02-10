from lib_plt import *
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

if len(sys.argv) < 2:
	print('Please supply filename')
	sys.exit()

filename = sys.argv[1]
data_location=os.path.join('data',filename)
confname = filename.split('.npz')[0]


def Data():
    global data, data1, data2, data3, data4, omega, ft_dos
    
    data = np.load(data_location, allow_pickle=True)
    ft_dos, omega = data['ft_dos'], data['omega']
    
    ft_dos=np.average(ft_dos,0)
    
    data1=ft_dos[:,:,0,0]
    data2=ft_dos[:,:,1,0] 
    data3=ft_dos[:,:,0,1] 
    data4=ft_dos[:,:,1,1]
    return    
Data()
##########################################################
######################## Main ############################
##########################################################
def Main():
    Data()
    data_length=np.shape(ft_dos)[0]
    n_x=np.shape(ft_dos)[0]
    n_y=np.shape(ft_dos)[1]
    n_spins=np.shape(ft_dos)[2]
    n_orbitals=np.shape(ft_dos)[3]
    
    text1 = (r'$\sigma=\uparrow$, orbital$=+$')    
    text2 = (r'$\sigma=\downarrow$, orbital$=+$')    
    text3 = (r'$\sigma=\uparrow$, orbital$=-$')    
    text4 = (r'$\sigma=\downarrow$, orbital$=-$')    
                   
    # max/min without central bright spot and lines:
    r=centre([0,0],n_x,n_y)
    r=(* r,)
    temp=ft_dos
    temp[r]=temp[0,0]
    temp[r[0],:]=temp[0,0]
    temp[:,r[1]]=temp[0,0]
    vmax=np.amax(temp)
    vmin=np.amin(temp)

    xlabel='$k_x a$'
    ylabel='$k_y a$'
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    # newcmp = ListedColormap(cm.afmhot(np.linspace(0, 0.7, 256)))
    newcmp = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black","orange","white"])
    
    eV=np.real(omega)
    epsilon=np.imag(omega)
    
    text_DOS = ('$|$FT(DOS)$|$'
    '\n'
    f'$\omega={eV:.0f}$'
    '\n'
    f'$\epsilon={epsilon}$'
    )    
    interpolation = 'none'
    #
    extent = [-n_x/2,n_x/2,-n_y/2,n_y/2]
    
    dos_map=data1
    # vmin, vmax = VminMax(data4)
    im=ax1.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)    
    ax1.text(0.67,0.05, 
        text1,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax1.transAxes)
    
    dos_map=data2
    ax2 = fig.add_subplot(222, sharey=ax1)
    ax2.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)
    ax2.text(0.67,0.05, 
        text2,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax2.transAxes)

    dos_map=data3
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
    # cbar.ax.set_yticklabels([vmin,vmax], color='white')
    cbar.ax.set_yticklabels(['Low','High'], color='white')
    ax3.text(0.67,0.05, 
        text3,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax3.transAxes)

    dos_map=data4
    ax4 = fig.add_subplot(224, sharex=ax2, sharey=ax1)
    im=ax4.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)
    ax4.text(0.67,0.05, 
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
    
    ax1.text(0.22,0.85, 
        text_DOS,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax1.transAxes)

    fig.set_size_inches(w=latex_width, h=latex_width) 

    plt.tight_layout()
    
    print('Done')
    return fig

latex_width=4.7747

fig=Main()
plt.show()
fig.savefig('out/'+confname+'.pdf', bbox_inches = "tight")