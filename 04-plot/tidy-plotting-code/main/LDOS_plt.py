from lib_plt import *

if len(sys.argv) < 2:
	print('Please supply name of configuration file')
	sys.exit()
    
data_location='..\\Data'

configuration = sys.argv[1]
confname = configuration.split('.conf')[0]
config_module = import_path(os.path.join(data_location,configuration))    
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

data = np.load(glob.glob(os.path.join(data_location,confname)+'.npz')[0], allow_pickle=True)
dos = data['dos']

dos=np.sum(dos,(2,3))
dos=dos[:,:,0,0,:]+dos[:,:,1,1,:]
##########################################################
########################## Main ##########################
##########################################################
def Main():
    index = FindNearestValueOfArray(np.real(omegas), omega)
    dos_map = dos[:,:,index]
    eV=np.real(omegas[index])
    epsilon=np.imag(omegas[index])
    
    text_DOS = ('DOS map'
    '\n'
    f'$\omega={eV:.2f}$'
    '\n'
    f'$\epsilon={epsilon}$'
    )    
    # k_F = Fermi_vector(mu, t)
    # lambda_F = Friedel_wavelength(k_F)
    # text_fermi = (f'$\lambda_F={lambda_F:.2f}$')
    # text = ('Model parameters: '
           # f'$\mu={mu:.2f}$, '
           # f'$s={s:.3f}$, '
           # f'$\delta={delta:.3f}$, '
           # f'$t={t:.2f}$, '
           # f'$d=({d[0]:.2f},{d[1]:.2f},{d[2]:.2f})$, '
           # f'$N_x = {n_sites_x}$, '
           # f'$N_y = {n_sites_y}$'
           # '\n'
           # 'Impurity location: '
           # f'{impurity_loc}, '
           # f'$V={V:.2f}$'
           # )
    interpolation = 'none'
    #
    vmin, vmax = np.amin(dos_map), np.amax(dos_map)
    dv=vmax-vmin
    vmax=vmin+dv#/50
    #
    extent = [-n_x/2,n_x/2,-n_y/2,n_y/2]
    #
    fig, axs = plt.subplots(1,1,figsize=(8,6.5))
    im = axs.imshow(
            dos_map.T, extent=extent,
            interpolation=interpolation,
            vmin=vmin, vmax=vmax,
            cmap=cm.YlOrRd)
    axs.set(xlabel='$x/a$', ylabel='$y/a$')
    
    fig.subplots_adjust(hspace=0.1, wspace=0.1)
    cbar = fig.colorbar(im, ax=axs)
#    cbar.ax.set_title(r'eV')
    fig.text(.167, .84, text_DOS,
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
    
omega=0
fig, axs = Main()
plt.show()
fig.savefig('out\\DOS_map_LaNiGa2.pdf', bbox_inches = "tight")