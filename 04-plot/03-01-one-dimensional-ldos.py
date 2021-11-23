from plt_lib import *
########################################################
########################## Data ########################
########################################################
energy=0
layer=(ALL,ALL)

omega=0
n_data=len(filenames)

ncol=1
nrow = n_data+1;

fig, axs = plt.subplots(nrows=nrow, ncols=ncol)

for i in range(n_data):
    filename = filenames[i]
    globals().update(conf_file(filename))

    import _pickle as cPickle
    with open(filename, 'rb') as f:
        data = cPickle.load(f)

    if phase=='lattice':

        j=0
        label_x = r'$x/a_x$'
        label_y = r'$y/a_y$'
        axs[0].set_xlabel(label_x)
        axs[0].set_ylabel(label_y)
        label_y=r'a)'

        plt_data = plot_data(fig,axs[j],data)

        fig,axs[j] = plt_data.lattice(energy=energy, model=data, layer=layer, spins=np.arange(n_spins), orbitals=np.arange(data.n_orbitals), annotate_orbs=False, show_cell_borders=True, show_basis_vectors=False, show_hopping=True, show_impurities=True, show_only_centre=True)

    else:
        if phase=='trivial':    
            j=2
            label_y=r'b)'
            axs[j].set_title(' ')
        if phase=='topological':
            j=3
            label_y=r'c)'
        if phase=='trivial_t':    
            j=4
            label_y=r'd)'
        if phase=='topological_t':
            j=5
            label_y=r'e)'
        if phase=='trivial_tt':
            j=6
            label_y=r'f)'
        if phase=='topological_tt':
            j=7
            label_y=r'g)'
        if phase=='trivial_ttt':
            j=8
            label_y=r'h)'
        if phase=='topological_ttt':
            j=9
            label_y=r'i)'

        axs[j].axis("off")

        plt_data = plot_data(fig,axs[j],data)

        fig, axs[j] = plt_data.differential_current_map(energy, layer=layer, orbital=[0,1], spin=None, cartesian=True, n_pts=40, gaussian_mean=0.2)

    text_box = AnchoredText(label_y, loc=2, pad=0.3, borderpad=0)
    plt.setp(text_box.patch, facecolor='white', alpha=0.3)
    axs[j].add_artist(text_box)
    fig.set_size_inches(w=latex_width, h=5) 
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
    label.append(rf'SSH model with additional nearest-neighbour intraorbital hopping $-t$. {label} This model shows that the topological state is delicate. Indeed, any nearest-neighbour intraorbital pairing model gives rise to a zero-energy state with density in the bulk.')
    label.append(r'a) Lattice cartoon with three sites, showing hopping terms, $v$, $w$ and $t$;')
    label.append(rf'local density of states at zero energy of an SSH model with {n_x} sites: '+
            r'b) Trivial phase: $v>w,t=0$;')
    label.append(r'c) Topological phase: $w>v,t=0$;')
    label.append(r'd) Trivial with hopping: $v>w>t \neq 0$;')
    label.append(r'e) Topological with hopping: $w>v>t \neq 0$;')
    label.append(r'f) Trivial with strong hopping: $v>t>w$;')
    label.append(r'g) Topological with strong hopping: $w>t>v$;')
    label.append(r'h) Trivial with stronger hopping: $t>v>w$;')
    label.append(r'i) Topological with stronger hopping: $t>w>v$.')
    text=' '.join(label)
    with open(output+'.txt', 'w') as f:
        f.write(text)

##########################################################
########################## Main ##########################
########################################################## 
main()
plt.savefig(output+'.pdf', bbox_inches = "tight")
caption()
