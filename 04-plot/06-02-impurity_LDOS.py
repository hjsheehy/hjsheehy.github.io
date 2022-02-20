from plt_lib import *

FILENAME='06-02-phase-diagrams.py'
FILENAMES=matrix_filenames()
print(FILENAMES)

def main():

    xlabel=r'x'
    ylabel=r'y'
    interpolation='none'

    cmap='plasma'
    (nx,ny)=np.shape(FILENAMES)
    nz=1
    data=np.zeros([nx,ny,nz])
    convergence=[]
    minmax_x=[]
    minmax_y=[]

    fig, axs = plt.subplots(1, 1, sharex='all', sharey='all')
    
    for filename in np.array(FILENAMES).flatten():
        globals().update(conf_file(filename))

        import _pickle as cPickle
        with open(filename, 'rb') as f:
            model = cPickle.load(f)

        magnetism=model.mean_magnetism(energy=energy)
        if magnetism>0.05:
            break

    # fig, axs = model.plot_lattice(fig, axs, atoms=None, plot_ldos=False, plot_magnetism=True, scale_factor=100)
    fig, axs = model.plot_lattice(fig, axs, atoms=None, plot_ldos=True, plot_magnetism=False, scale_factor=10)

    i=j=k=0
    # im = axs[i,j].imshow(data[k].T, interpolation=interpolation, cmap=cmap[k], origin='lower', extent=extent)
    # axs[i,j].set(title=r'LDOS')
    # divider = make_axes_locatable(axs[i,j])
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    
    # plt.colorbar(im, cax=cax)

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH/10) 

def caption():
    label=[]
    label.append(r'The phase diagrams for an SSH model with multiorbital, equal-spin superconductivity interaction and no impurity. Two initial pairing fields are given to the self-consistent loop: non-unitary triplet and equal-spin triplet. The solution with the lowest free energy is selected. We observe that the mean-field spin-sectors are equal $\overline{\Delta_{v/w}^{\uparrow\uparrow}}=\overline{\Delta_{v/w}^{\downarrow\downarrow}}$, hence the equal-spin triplet theory is favoured.')
    text=' '.join(label)
    with open(title+'.txt', 'w') as f:
        f.write(text)

title=OUTPUT+'-no_impurity'
# applied bias for e.g. ldos, staggered density and IPR
energy=0

# caption()
main()
# plt.savefig(title+'.pdf', bbox_inches = "tight")
plt.show()
