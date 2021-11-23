from plt_lib import *
########################################################
########################## Data ########################
########################################################
energy=0
layer=(ALL,ALL,0)

omega=0
n_data=len(filenames)

nrow = 1; ncol = n_data;

fig, axs = plt.subplots(nrows=nrow, ncols=ncol)

for i in range(n_data):
    filename = filenames[i]
    globals().update(conf_file(filename))

    import _pickle as cPickle
    with open(filename, 'rb') as f:
        data = cPickle.load(f)

    plot_data = plot_data(fig,axs[i],data)

    fig,axs[i] = plot_data.lattice(energy=energy, model=data, layer=layer, spins=np.arange(n_spins), orbitals=np.arange(data.n_orbitals), annotate_orbs=False, show_cell_borders=True, show_basis_vectors=False, show_hopping=True, show_impurities=True, show_only_centre=True)

########################################################
########################## Plot ########################
########################################################
def main():
    label_x = r'$x/a_x$'
    label_y = r'$y/a_y$'
    
    i=0
    axs[i].set_ylabel(label_y)
    i=1
    plt.setp(axs[i].get_yticklabels(), visible=False)

    for i in range(n_data):
        axs[i].set_xlabel(label_x)

    plt.tight_layout()

    fig.set_size_inches(w=latex_width, h=4.5) 
    return fig
########################################################
######################### Caption ######################
########################################################
def caption():
    fermi_vector = Fermi_vector(mu, t, omega)
    friedel_wavelength = Friedel_wavelength(fermi_vector)
    
    text=rf'''Two simple lattices with ${n}\times{n}$ lattice sites, with single (blue) orbital (right) or multiorbital (red, blue) (left). 
Lattice impurities are gives by the black dot with a central asterix. 
Hopping vectors for the first site are denoted by black, dashed arrows.
'''
    with open(output+'.txt', 'w') as f:
        f.write(text)

##########################################################
########################## Main ##########################
########################################################## 
main()
plt.savefig(output+'.pdf', bbox_inches = "tight")
caption()
