from plt_lib import *

FILENAMES=matrix_filenames()

def main():

    xlabel=r'$U_v$'
    ylabel=r'$U_w$'
    interpolation='none'

    data_dict={
            'Free energy': 0,
            'Magnetism': 1,
            r'$\overline{|\Delta_v^{\uparrow\uparrow}|}$': 2,
            r'$\overline{|\Delta_w^{\uparrow\uparrow}|}$': 3,
            r'$\overline{|\Delta_v^{\downarrow\downarrow}|}$': 4,
            r'$\overline{|\Delta_w^{\downarrow\downarrow}|}$': 5,
            r'Staggered density': 6,
            r'IPR': 7,
            }
    cmap=['plasma',
            'bwr',
            'gnuplot',
            'Greens',
            'gnuplot',
            'Greens',
            'Reds',
            'Blues']
    (nx,ny)=np.shape(FILENAMES)
    nz=len(data_dict)
    data=np.zeros([nx,ny,nz])
    convergence=[]
    minmax_x=[]
    minmax_y=[]

    fig, axs = plt.subplots(4, 2, sharex='all', sharey='all')
    
    for i,row in enumerate(FILENAMES):
        for j,filename in enumerate(row):
            globals().update(conf_file(filename))

            import _pickle as cPickle
            with open(filename, 'rb') as f:
                model = cPickle.load(f)

            minmax_x.append(U_v)
            minmax_y.append(U_w)

            data[i,j,0]=model.free_energy[-1]
            data[i,j,1]=model.mean_magnetism(energy=energy)
            data[i,j,2]=np.mean(model.gorkov(atom_i='A', atom_f='B', orbital_i=None, orbital_f=None, hop_vector=None, spin_i='up', spin_f='up'))
            data[i,j,3]=np.mean(model.gorkov(atom_i='B', atom_f='A', orbital_i=None, orbital_f=None, hop_vector=[1,0], spin_i='up', spin_f='up'))
            data[i,j,4]=np.mean(model.gorkov(atom_i='A', atom_f='B', orbital_i=None, orbital_f=None, hop_vector=None, spin_i='dn', spin_f='dn'))
            data[i,j,5]=np.mean(model.gorkov(atom_i='B', atom_f='A', orbital_i=None, orbital_f=None, hop_vector=[1,0], spin_i='dn', spin_f='dn'))
            data[i,j,6]=model.mean_staggered_density(atom_i='A',atom_f='B',energy=energy)
            data[i,j,7]=model.IPR_staggered_density(atom_i='A',atom_f='B',energy=energy)
            # data[i,j,5]=model.converged
            if not model.converged:
                convergence.append([U_v,U_w])

    min_x=min(minmax_x)
    max_x=max(minmax_x)
    min_y=min(minmax_y)
    max_y=max(minmax_y)
    extent=[min_x,max_x,min_y,max_y]

    k=0
    for i in range(4):
        for j in range(2):
            im = axs[i,j].imshow(data[...,k].T, interpolation=interpolation, cmap=cmap[k], origin='lower', extent=extent)
            axs[i,j].set(title=list(data_dict)[k])
            divider = make_axes_locatable(axs[i,j])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            
            if convergence==[]:
                x,y=[],[]
            else:
                x,y=list(np.transpose(convergence))
            axs[i,j].scatter(x=x, y=y, c='w', marker='s', s=10)
            axs[i,j].set_xlim([min_x,max_x])
            axs[i,j].set_ylim([min_y,max_y])
            plt.colorbar(im, cax=cax)

            k+=1

    # axes labels:
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    fig.set_size_inches(w=LATEX_WIDTH, h=6) 
    plt.savefig(title+'.pdf', bbox_inches = "tight")

def caption():
    label=[]
    label.append(r'The phase diagrams for an SSH model with multiorbital, equal-spin superconductivity interaction and no impurity. Two initial pairing fields are given to the self-consistent loop: non-unitary triplet and equal-spin triplet. The solution with the lowest free energy is selected. We observe that the mean-field spin-sectors are equal $\overline{\Delta_{v/w}^{\uparrow\uparrow}}=\overline{\Delta_{v/w}^{\downarrow\downarrow}}$, hence the equal-spin triplet theory is favoured.')
    text=' '.join(label)
    with open(title+'.txt', 'w') as f:
        f.write(text)

title=OUTPUT+'-no_impurity'
# applied bias for e.g. ldos, staggered density and IPR
energy=0

caption()
main()
