from plt_lib import *
########################################################
########################## Data ########################
########################################################
i=0
filename = filenames[i]
globals().update(conf_file(filename))
import _pickle as cPickle
with open(filename, 'rb') as f:
    data = cPickle.load(f)
########################################################
########################## Plot ########################
########################################################
import matplotlib.pyplot as plt
latex_width=4.7747
n_pts=40    

# fig,axs = plot_data.lattice(energy=energy, model=data, layer=layer, spins=np.arange(plot_data.n_spins), orbitals=np.arange(data.n_orbitals), annotate_orbs=True, show_cell_borders=True, show_basis_vectors=False, show_hopping=True, show_impurities=True, show_only_centre=True)

def fig_a():
    energy=0

    fig, axs = plt.subplots(1)
    plt_data = plot_data(fig,axs,data)

    fig, axs = plt_data.differential_current_map(energy, layer=(ALL,ALL), orbital=[0,1], spin=None, cartesian=True, n_pts=n_pts, gaussian_mean=0.1)

    label_x = r'$x/a_x$'
    label_y = r'$y/a_y$'
    axs.set_ylabel(label_y)
    axs.set_xlabel(label_x)

    plt.tight_layout()

    axs.set_title('zero-bias local density of states')
    return fig

def fig_b():
    dimension=0

    fig, axs = plt.subplots(1)
    plt_data = plot_data(fig,axs,data)

    fig,axs = plt_data.band_structure(dimension)
    axs.set_title('Bandstructure')
    fig.set_size_inches(w=latex_width, h=4.5) 
    return fig

def fig_c():
    energy=0
    colors = [cm.rainbow(i) for i in np.linspace(0, 1, 3)]

    fig, axs = plt.subplots(1)

    centre=data.centre
    ldos_a = data.local_density_of_states(energy, orbital=0)
    ldos_b = data.local_density_of_states(energy, orbital=1)
    ldos=ldos_a-ldos_b
    ldos=ldos[:,centre[1]]
    vmin,vmax=np.amin(ldos),np.amax(ldos)
    gorkov=data.gorkov()

    y=centre[1]
    Delta_T=np.real([[gorkov[x,y,0,s,x,y,1,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
    Delta_T=np.sum(Delta_T,0)/data.n_spins
    Delta_R=np.real([[gorkov[x-1,y,1,s,x,y,0,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
    Delta_R=np.sum(Delta_R,0)/data.n_spins

    x=np.arange(-centre[0],centre[0]+1)
    for k,y in enumerate([ldos,Delta_T,Delta_R]):
        y=np.fft.fftshift(y)
        axs.plot(x,y, 
                color=colors[k])
    axs.legend([r'CDW', r'$\Delta_T$',r'$\Delta_R$'],labelcolor='k')

    label_x = r'$x/a_x$'
    label_y = r'Field strength'
    axs.set_ylabel(label_y)
    axs.set_xlabel(label_x)

    plt.tight_layout()

    axs.set_title('CDW, INT and repulsion')
    fig.set_size_inches(w=latex_width, h=4.5) 
    return

fig_a()
plt.savefig(output+'_a.pdf', bbox_inches = "tight")
plt.close()

fig_b()
plt.savefig(output+'_b.pdf', bbox_inches = "tight")
plt.close()

fig_c()
plt.savefig(output+'_c.pdf', bbox_inches = "tight")
