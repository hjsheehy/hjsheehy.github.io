from plt_lib import *
########################################################
########################## Data ########################
########################################################
for filename in filenames:
    globals().update(conf_file(filename))
    import _pickle as cPickle
    with open(filename, 'rb') as f:
        data = cPickle.load(f)
    ########################################################
    ########################## Plot ########################
    ########################################################
    import matplotlib.pyplot as plt
    n_pts=40    

    # fig,axs = plot_data.lattice(energy=energy, model=data, layer=layer, spins=np.arange(plot_data.n_spins), orbitals=np.arange(data.n_orbitals), annotate_orbs=True, show_cell_borders=True, show_basis_vectors=False, show_hopping=True, show_impurities=True, show_only_centre=True)

    fig, axs = plt.subplots(3,figsize=(LATEX_WIDTH, 7.5), gridspec_kw={'width_ratios':[1], 'height_ratios':[20,20,1]})

    energy=0
    dimension=0

    plt_data = plot_data(fig,axs[2],data)

    fig, axs[2] = plt_data.differential_current_map(energy, layer=(ALL,ALL), orbital=[0,1], spin=None, cartesian=True, n_pts=n_pts, gaussian_mean=0.1)

    label_x = r'$x/a_x$'
    label_y = r'$y/a_y$'
    axs[2].set_ylabel(label_y)
    axs[2].set_xlabel(label_x)

    axs[2].set_title('Zero-bias local density of states')

    plt_data = plot_data(fig,axs[1],data)

    fig,axs[1] = plt_data.band_structure(dimension)
    axs[1].set_title('Bandstructure')
    label_x = r''
    axs[1].set_xlabel(label_x)

    colors = [cm.rainbow(i) for i in np.linspace(0, 1, 4)]

    centre=data.centre
    ldos_a = data.local_density_of_states(energy, orbital=0)
    ldos_b = data.local_density_of_states(energy, orbital=1)
    stag=ldos_a-ldos_b
    ldos=ldos_a+ldos_b
    stag=stag[:,centre[1]]
    ldos=ldos[:,centre[1]]
    vmin,vmax=np.amin(ldos),np.amax(ldos)
    gorkov=data.gorkov()

    y=centre[1]
    Delta_T=np.real([[gorkov[x,y,0,s,x,y,1,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
    Delta_T=np.sum(Delta_T,0)/data.n_spins
    Delta_R=np.real([[gorkov[x-1,y,1,s,x,y,0,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
    Delta_R=np.sum(Delta_R,0)/data.n_spins

    x=np.arange(-centre[0],centre[0]+1)
    for k,y in enumerate([ldos,stag,Delta_T,Delta_R]):
        y=np.fft.fftshift(y)
        axs[0].plot(x,y, 
                color=colors[k])
    axs[0].legend([r'Zero-bias LDOS',r'Staggered density', r'$\Delta_T$',r'$\Delta_R$'],labelcolor='k',loc='upper right')

    label_x = r''
    label_y = r'Field strength'
    axs[0].set_ylabel(label_y)
    axs[0].set_xlabel(label_x)

    axs[0].set_title(f'$U_T={U_T}$, $U_R={U_R}$, $\mu={mu}, V={V}$, {state}')

    plt.tight_layout()

    plt.savefig(output+f'_{U_T}_{U_R}.pdf', bbox_inches = "tight")
