
from plt_lib import *
########################################################
########################## Data ########################
########################################################
energy=0
layer=(ALL,ALL)

omega=0
n_data=len(filenames)

ncol=1
nrow = n_data

fig, axs = plt.subplots(nrows=nrow, ncols=ncol)

for i in range(n_data):
    filename = filenames[i]
    globals().update(conf_file(filename))

    import _pickle as cPickle
    with open(filename, 'rb') as f:
        data = cPickle.load(f)

    # if phase=='lattice':
    #     j=0
    #     label_x = r'$x/a_x$'
    #     label_y = r'$y/a_y$'
    #     axs[0].set_xlabel(label_x)
    #     axs[0].set_ylabel(label_y)
    #     label_y=r'a)'

    #     plt_data = plot_data(fig,axs[j],data)

    #     fig,axs[j] = plt_data.lattice(energy=energy, model=data, layer=layer, spins=np.arange(n_spins), orbitals=np.arange(data.n_orbitals), annotate_orbs=False, show_cell_borders=True, show_basis_vectors=False, show_hopping=True, show_impurities=True, show_only_centre=True)

#    else:
    if phase=='spin_chain':
        j=0
        label_y=r'b)'
        axs[j].set_title(' ')
    if phase=='small_B':    
        j=1
        label_y=r'c)'
    if phase=='large_B':
        j=2
        label_y=r'd)'
    if phase=='small_SC':
        j=3
        label_y=r'e)'
    if phase=='large_SC':
        j=4
        label_y=r'f)'
    if phase=='small_SC_small_B':
        j=5
        label_y=r'g)'
    if phase=='large_SC_small_B':
        j=6
        label_y=r'h)'
    if phase=='small_SC_large_B':
        j=7
        label_y=r'i)'
    if phase=='large_SC_large_B':
        j=8
        label_y=r'j)'

    axs[j].axis("off")

    plt_data = plot_data(fig,axs[j],data)

    fig, axs[j] = plt_data.differential_current_map(energy, layer=layer, orbital=None, spin=None, cartesian=True, n_pts=40, gaussian_mean=0.2)

    text_box = AnchoredText(label_y, loc=3, pad=0.3, borderpad=0)
    plt.setp(text_box.patch, facecolor='white', alpha=0.3)
    axs[j].add_artist(text_box)
    fig.set_size_inches(w=LATEX_WIDTH, h=5) 
    axs[1].axis("off")
########################################################
########################## Plot ########################
########################################################
def main():
    
    #i=0
    #i=1
    #plt.setp(axs[i].get_yticklabels(), visible=False)

    #for i in range(n_data):
    #    axs[i].set_xlabel(label_x)

    return fig
########################################################
######################### Caption ######################
########################################################
def caption():
    label=[]
    label.append(rf'Tight-binding chain with two orbitals, $A$ (blue), $B$ (red), with chemical potentials $mu$ and $-mu$, and hopping parameter $-t$. The model exhibits two qualitatively distinct topological phases, depending on the value of $mu$. We refer to these phases as trivial and topological in reference to the SSH model, although we observe that localised modes are ubiquitous, albeit, bulk density is lowest in the topological phase.')
    label.append(rf'a) Lattice cartoon with three sites, showing the hopping term a$t$={{t}};')
    label.append(rf'local density of states at zero energy of an SSH model with {n_x} sites: '+
            rf'b) Trivial phase: $\mu={mu}$;')
    label.append(rf'c) Topological phase: $\mu={mu}$;')
    label.append(rf'd) Trivial phase: $\mu={mu}$;')
    label.append(rf'e) Topological phase: with $\mu={mu}$;')
    label.append(rf'f) Trivial phase: $\mu={mu}$;')
    label.append(rf'g) Topological phase: $\mu={mu}$;')
    label.append(rf'h) Trivial phase: $\mu={mu}$;')
    text=' '.join(label)
    with open(output+'.txt', 'w') as f:
        f.write(text)

##########################################################
########################## Main ##########################
########################################################## 
main()
plt.savefig(output+'.pdf', bbox_inches = "tight")
caption()
