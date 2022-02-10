from plt_lib import *

FILENAMES=matrix_filenames()
xlabel=r'$U_v$'
ylabel=r'$U_w$'
interpolation='none'

def plot_phase_diagram():
    energy=0
    data_dict={
            'Free energy': 0,
            'Magnetism': 1,
            r'$\overline{|\Delta_T^{\uparrow\uparrow}|}$': 2,
            r'$\overline{|\Delta_R^{\uparrow\uparrow}|}$': 3,
            r'Staggered density': 4,
            r'IPR': 5,
            }
    cmap=['plasma',
            'bwr',
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

    fig, axs = plt.subplots(3, 2, sharex='all', sharey='all')
    
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
            data[i,j,2]=np.mean(model.gorkov(atom_i='A', atom_f='B', orbital_i=None, orbital_f=None, hop_vector=None, spin_i=None, spin_f=None))
            data[i,j,3]=np.mean(model.gorkov(atom_i='B', atom_f='A', orbital_i=None, orbital_f=None, hop_vector=[1,0], spin_i='down', spin_f='up'))
            data[i,j,4]=model.mean_staggered_density(atom_i='A',atom_f='B',energy=energy)
            data[i,j,5]=model.IPR_staggered_density(atom_i='A',atom_f='B',energy=energy)
            # data[i,j,5]=model.converged
            if not model.converged:
                convergence.append([U_v,U_w])
    min_x=min(minmax_x)
    max_x=max(minmax_x)
    min_y=min(minmax_y)
    max_y=max(minmax_y)
    extent=[min_x,max_x,min_y,max_y]

    k=0
    for i in range(3):
        for j in range(2):
            im = axs[i,j].imshow(data[...,k].T, interpolation=interpolation, cmap=cmap[k], origin='lower', extent=extent)
            if i==0:
                axs[i,j].set(xlabel=xlabel)
            if j==2:
                axs[i,j].set(ylabel=ylabel)
            axs[i,j].set(title=list(data_dict)[k])
            divider = make_axes_locatable(axs[i,j])
            cax = divider.append_axes("right", size="5%", pad=0.05)

            x,y=list(np.transpose(convergence))
            axs[i,j].scatter(x=x, y=y, c='w', marker='s', s=10)
            axs[i,j].set_xlim([min_x,max_x])
            axs[i,j].set_ylim([min_y,max_y])
            plt.colorbar(im, cax=cax)

            k+=1


    fig.set_size_inches(w=LATEX_WIDTH, h=5) 
    plt.savefig(OUTPUT+'SSH-triplet_with_impurity.pdf', bbox_inches = "tight")

def plot_perculation():
    energy=0
    data_dict={
            'Seed': 0,
            'Layer': 1,
            }

    (nx,ny)=np.shape(FILENAMES)
    nz=len(data_dict)
    data=np.zeros([nx,ny,nz])
    for i,row in enumerate(FILENAMES):
        for j,filename in enumerate(row):
            globals().update(conf_file(filename))

            import _pickle as cPickle
            with open(filename, 'rb') as f:
                model = cPickle.load(f)

            data[i,j,0]=filename_seed(filename)
            data[i,j,1]=filename_layer(filename)

    for i in range(2):
        fig = plt.figure()
        axs = fig.add_subplot(111)
        im = axs.imshow(data[...,i],interpolation=interpolation)
        axs.set(xlabel=xlabel, ylabel=ylabel)
        # axs.set(title=list(data_dict)[i])
        values = np.unique(data[...,i].ravel()).astype(int)
        colors = [ im.cmap(im.norm(value)) for value in values]
        # create a patch (proxy artist) for every color 
        patches = [ mpatches.Patch(color=colors[j], label=f"{list(data_dict)[i]} {values[j]}") for j in range(len(values)) ]
        # put those patched as legend-handles into the legend
        plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
        plt.show()
        plt.close()

plot_phase_diagram()

# plot_perculation()
