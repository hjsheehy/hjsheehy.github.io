from lib_plt import *
import matplotlib.ticker as mtick

if len(sys.argv) < 2:
	print('Please supply filename')
	sys.exit()

# data_location=os.path.join('data','superconducting_state')
data_location='data'

configuration = sys.argv[1]
confname = configuration.split('.conf')[0]
config_module = import_path(os.path.join(data_location,configuration))    
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

data = np.load(os.path.join(data_location,'DOS_'+confname+'.npz'), allow_pickle=True)
dos = data['dos']
##########################################################
######################## Main ############################
##########################################################
def Main():
    index = FindNearestValueOfArray(np.real(omegas), omega)
    dos_map = dos[:,:,:,:,index]
    dos_map=dos_map.sum((2,3))
    dos_ft=FT(dos_map)

    eV=np.real(omegas[index])
    epsilon=np.imag(omegas[index])

    text_DOS = ('FT(DOS)'
    '\n'
    f'$\omega={eV:.2f}$'
    '\n'
    f'$\epsilon={epsilon}$'
    )    
    interpolation = 'none'
    # max/min without central bright spot and lines:
    temp=dos_map
    r=centre([0,0],n_x,n_y)
    r=(* r,)
    temp[r]=temp[0,0]
    # temp[r[0],:]=temp[0,0]
    # temp[:,r[1]]=temp[0,0]
    vmax=np.amax(temp)
    vmin=np.amin(temp)
    # dv=vmax-vmin
    # vmax=vmin+dv
    
    # dos_map[r]=vmax
    #
    extent = [-n_x/2,n_x/2,-n_y/2,n_y/2]
    #
    fig, axs = plt.subplots(1,1,figsize=(latex_width,6.5))
    # newcmp = ListedColormap(cm.gist_earth(np.linspace(0, 0.75, 256)))
    newcmp = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black","orange","white"])
    im = axs.imshow(
            dos_map.T, extent=extent,
            interpolation=interpolation,
            vmin=vmin, vmax=vmax,
            cmap=newcmp)
    axs.set(xlabel='$k_x a$', ylabel='$k_y a$')
    axs.set_xticks([-n_x/2,0,n_x/2])
    axs.set_xticklabels(['$-\pi$',0,'$\pi$'])
    axs.set_yticks([-n_y/2,0,n_y/2])
    axs.set_yticklabels(['$-\pi$',0,'$\pi$'])
    
    fig.subplots_adjust(hspace=0.1, wspace=0.1)
    cbar = fig.colorbar(im, ax=axs, ticks=[vmin,vmax])
    cbar.ax.set_yticklabels(['Low','High']) 
#    cbar.ax.set_title(r'eV')
    fig.text(.17, .835, text_DOS,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center')
    # fig.text(0.7, .86, text_fermi,
             # {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           # ec="none", pad=0.2)}, ha='center', va='center')

    # fig.text(.45, -0.03, text,
             # {'bbox': dict(boxstyle="round", fc="white", ec="black", pad=0.3)},
             # ha='center', va='bottom')
    #
    return fig, axs

latex_width=4.7747
omega=0

fig, ax = Main()
plt.show()
fig.savefig('out/'+'LDOS_FT_'+confname+'.pdf', bbox_inches = "tight")