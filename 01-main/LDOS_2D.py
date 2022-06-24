from lib import *

FILENAME=sys.argv[0].split('.')[0]
DATA=os.path.join(DATA,FILENAME)
FIG=os.path.join(FIG,FILENAME)
for directory in [FIG,DATA]:
    if not os.path.exists(directory):
        os.makedirs(directory)

def main(CALCULATE=False):

    DATA_NAME='single_impurity'
    DATA_NAME=os.path.join(DATA,DATA_NAME)

    A=Atom([0,0],'A')
    A.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    tb=TightBinding(lattice_vectors,'SSH')
    tb.add_atom(A)
    tb.n_spins=1

    tb.cut(n_cells, [0,1], glue_edgs=True)
    tb.set_onsite(-mu,orbital='s')
    tb.set_hopping(-t,hop_vector=[1,0],label='t')
    tb.set_hopping(-t,hop_vector=[0,1],label='t')
    
    tb.add_impurities(V,impurity_locations,label='V')

    if CALCULATE:
        tb.solve()

        energy_interval=np.linspace(0,0,1)
        resolution=0.1
        greens_function_xy=GreensFunction(tb,energy_interval,resolution, k_axes=None)

        with open(DATA_NAME, 'wb') as f:
            cPickle.dump(greens_function_xy, f)
    else:
        greens_function_xy = np.load(DATA_NAME, allow_pickle=True)

    return greens_function_xy

def plot():

    FIGNAME='LDOS_single_impurity'

    cmap=cm.YlOrRd

    offset=10**(-6)
        
    Z=np.fft.fftshift(ldos)
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
    fig.set_size_inches(w=LATEX_WIDTH, h=0.9*LATEX_WIDTH) 
    ax.set(title=title,
    xlabel=r"$r_x/a$", ylabel=r"$r_y/a$")
    # Text label bug: hidden by contour
    # ax.text(-1, 1, -7, 'FS', color='black')
    # tick labels
    n = greens_function_xy._pieces[0]
    c=int(n/2)
    label_loc=[0,c,n]
    labels=[rf'$-{c}$',r'$0$',rf'${c}$']
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

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    plt.tight_layout()
    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight", dpi=DPI)
    return fig, ax


############################################################

n_cells=43
mu=-3.7
t=1
V=.1
impurity_locations=[[0,0]]
impurity_locations=points_on_circle(radius=9, x0=0, y0=0)

greens_function_xy = main(CALCULATE=True)

energy=0
ldos = greens_function_xy.local_density_of_states(energy=energy)

elev=35
azim=135
title=r'$-\frac{1}{\pi}\Im\hat{G}^R(\omega+i\epsilon;\mathbf{r},\mathbf{r})$'+r'$|_{{\omega={}}}$'.format(energy)

fig,ax = plot()

plt.show()
