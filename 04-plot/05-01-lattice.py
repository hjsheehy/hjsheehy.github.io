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
fig, axs = plt.subplots(1)
latex_width=4.7747
fig.set_size_inches(w=latex_width, h=4.5) 

energy=0
layer=(slice(None),slice(None))
plot_data = plot_data(fig,axs,data)
#dos = plot_data.density_of_states
orbital=0

fig,axs = plot_data.lattice(energy=energy, model=data, layer=layer, spins=np.arange(plot_data.n_spins), orbitals=np.arange(data.n_orbitals), annotate_orbs=True, show_cell_borders=True, show_basis_vectors=False, show_hopping=True, show_impurities=True, show_only_centre=True)


def main():
    label_x = r'$x/a_x$'
    label_y = r'$y/a_y$'
    
    i=0
    axs.set_ylabel(label_y)

    axs.set_xlabel(label_x)

    plt.tight_layout()
    return fig
#fig = main()
#energy
dimension=0
#fig,axs = plot_data.band_structure(dimension)

#fig, axs = plot_data.differential_current_map(energy, layer=(ALL,ALL), orbital=[0,1], spin=None, cartesian=True, n_pts=40, gaussian_mean=0.2)
#plt.show()
plt.savefig(output+'.pdf', bbox_inches = "tight")
