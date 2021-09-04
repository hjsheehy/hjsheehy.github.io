from plt_lib import *
########################################################
########################## Data ########################
########################################################
layer=(slice(None),slice(None),0)
omega=0
n_data=len(filenames)

nrow = 1; ncol = n_data;

fig, axs = plt.subplots(nrows=nrow, ncols=ncol)

for i in range(n_data):
    filename = filenames[i]
    globals().update(conf_file(filename))

    import _pickle as cPickle
    with open(filename, 'rb') as f:
        model = cPickle.load(f)
    ldos = LDOS(model.density_of_states, model.omegas, omega, trace_over=False, layer=layer)
    lattice = Lattice(fig, axs[i], model=model, ldos=ldos, layer=layer, spins=np.arange(n_spins), orbs=np.arange(n_orbs), annotate_orbs=False, show_cell_borders=True, show_basis_vectors=False, show_hopping=True, show_impurities=True, show_only_centre=True)

    axs[i] = lattice.ax
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

    text=rf'''
Two simple lattices with ${n}\times{n}$ lattice sites, with single (blue) orbital (right) or multiorbital (red, bleu) (left). 
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
