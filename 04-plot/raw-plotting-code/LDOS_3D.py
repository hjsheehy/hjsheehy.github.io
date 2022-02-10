from lib_plt import *

if len(sys.argv) < 2:
	print('Please supply name of configuration file')
	sys.exit()
    
data_location='data'

data_location=os.path.join(data_location,'normal_state')

configuration = sys.argv[1]
confname = configuration.split('.conf')[0]
config_module = import_path(os.path.join(data_location,configuration))    
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

data = np.load(glob.glob(os.path.join(data_location,'DOS_'+confname)+'.npz')[0], allow_pickle=True)
dos = data['dos']

dos = dos.sum((2,3))
###############################################
def Main():
    offset=10**(-6)
    #######################################
    ################# Data ################
    #######################################
    index = FindNearestValueOfArray(np.real(omegas), omega)

    LDOS = dos[:,:,index]

    Z=LDOS
    X,Y = np.mgrid[:Z.shape[0],:Z.shape[1]]
    #######################################
    ############# Plot data ###############
    #######################################
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    #3D:
    vmax=np.amax(Z)
    vmin=np.amin(Z)
    dv=vmax-vmin
    vmin=vmax-dv/3
    ax.plot_surface(X, Y, Z, cmap=cmap, lw=0, rstride=1, cstride=1, linewidth=0,antialiased=False,
        vmax=vmax, vmin=vmin)
    #2D projection:
    lower_2D_projection = 0.5 #8
    min=vmin-vmin*lower_2D_projection
    Z=offset*Z+min
    vmax=np.amax(Z)
    vmin=np.amin(Z)
    dv=vmax-vmin
    vmin=vmax-dv/3
    ax.contourf(X, Y, Z, 50, cmap=cmap, antialiased=True,
        vmax=vmax, vmin=vmin)
    # ax.contour( X, Y, offset*Z+min, 10, colors="k", linewidths=0.2, linestyles="solid")
    #######################################
    ######## Plot customisation ###########
    #######################################
    #width w from latex command \printinunitsof{in}\prntlen{\textwidth}
    fig.set_size_inches(w=latex_width, h=0.9*latex_width) 
    ax.set(title=title,
    xlabel=r"$r_x/a$", ylabel=r"$r_y/a$")
    # Text label bug: hidden by contour
    # ax.text(-1, 1, -7, 'FS', color='black')
    # tick labels
    label_loc=[0,21,43]
    labels=[r'$-21$',r'$0$',r'$21$']
    plt.xticks(label_loc, labels)
    plt.yticks(label_loc, labels)
    # ax.set(zlabel=r'LDOS')
    # ax.text(X.min()*1.1, Y.min()*1.2, Z.max()*1.3, r'$\epsilon_0(\mathbf{k})/t$')
    # ax.set_zlim(vmin, vmax)
    # Prevent label rotation
    ax.xaxis.set_rotate_label(False)
    ax.yaxis.set_rotate_label(False)
    ax.zaxis.set_rotate_label(False)
    # Clear grid & create an empty cubic wireframe 
    ax.grid(False)
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    # remove spacing between tick labels and axes3d
    [t.set_va('center') for t in ax.get_yticklabels()]
    [t.set_ha('left') for t in ax.get_yticklabels()]
    [t.set_va('center') for t in ax.get_xticklabels()]
    [t.set_ha('right') for t in ax.get_xticklabels()]
    [t.set_va('center') for t in ax.get_zticklabels()]
    [t.set_ha('left') for t in ax.get_zticklabels()]
    # ticks pointing into the plot that leave a clean edge around the boundary
    ax.xaxis._axinfo['tick']['inward_factor'] = 0
    ax.xaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.yaxis._axinfo['tick']['inward_factor'] = 0
    ax.yaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.zaxis._axinfo['tick']['inward_factor'] = 0
    ax.zaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.zaxis._axinfo['tick']['outward_factor'] = 0.4
    #
    # ax.xaxis.set_major_locator(MultipleLocator(1))
    # ax.yaxis.set_major_locator(MultipleLocator(1))
    # ax.zaxis.set_major_locator(MultipleLocator(4))
    ax.view_init(elev=elev, azim=azim)
    # fig.tight_layout()
    return fig, ax
###############################################
#################### Main #####################
###############################################
omega=0
latex_width=4.7747
elev=45
azim=135
title=r'$-\frac{1}{\pi}\Im\hat{G}^R(\omega+i\epsilon;\mathbf{r},\mathbf{r})$'+r'$|_{{\omega={}}}$'.format(omega)
# title=r'Local density of states'

# cmap='viridis'
# cmap='gist_heat'
cmap=cm.YlOrRd
fig,ax=Main()
plt.show()
fig.savefig('out/'+'LDOS_3D_'+confname + '.pdf')