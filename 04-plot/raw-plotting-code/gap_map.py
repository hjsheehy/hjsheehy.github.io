from lib_plt import *

if len(sys.argv) < 2:
	print('Please supply name of configuration file')
	sys.exit()
    
# data_location=os.path.join('data','spectrum_singlet')
data_location='..\\Data'

configuration = sys.argv[1]
confname = configuration.split('.conf')[0]
config_module = import_path(os.path.join(data_location,configuration))    
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

data = np.load(os.path.join(data_location,confname+'.npz'), allow_pickle=True)
dos = data['dos']

##########################################################
########################## Main ##########################
##########################################################
print('''Code adapted extensively from `Gapmap Algorithm' of Kristine Lang and subsequently extensively revised by Jennifer Eve Hoffman (see her Appendix A of her thesis)''')
def temp():
    interpolation = 'none'
    #
    # vmin, vmax = np.amin(dos_map), np.amax(dos_map)
    # dv=vmax-vmin
    # vmax=vmin+dv#/50
    #
    extent = [-n_x/2,n_x/2,-n_y/2,n_y/2]
    #
    fig, axs = plt.subplots(1,1,figsize=(8,6.5))
    im = axs.imshow(
            dos_map.T, extent=extent,
            interpolation=interpolation,
            # vmin=vmin, vmax=vmax,
            cmap=cm.YlOrRd)
    axs.set(xlabel='$x/a$', ylabel='$y/a$')
    
    fig.subplots_adjust(hspace=0.1, wspace=0.1)
    cbar = fig.colorbar(im, ax=axs)
    return fig, axs

en1, en2 = 0.5,1.5
skipresfil = True
skiprespct = 90
outrangefil = False
outrangepct = 80
contiguity = False
T1 = 0.1
badPixel = False

gapMap=GapMap(dos, omegas, en1, en2,
        skipresfil, skiprespct, outrangefil, outrangepct, contiguity, T1, badPixel,
        n_x, n_y, n_spins, n_orbitals)
dos_map=gapMap.sum((2,3))
fig, axs = temp()
plt.show()

fig.savefig('out\\gap_map_LaNiGa2_singlet.pdf', bbox_inches = "tight")