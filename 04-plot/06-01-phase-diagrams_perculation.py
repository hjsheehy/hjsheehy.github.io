from plt_lib import *

FILENAMES=matrix_filenames()
xlabel=r'$U_v$'
ylabel=r'$U_w$'
interpolation='none'
cmap=['prism','magma']

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

    xlabel=r'$U_v$'
    ylabel=r'$U_w$'
    title=OUTPUT+'-phase-diagram-perculation'
    fig, axs = plt.subplots(2, 1, sharex='all', sharey='all')
    for i in range(2):
        im = axs[i].imshow(data[...,i],interpolation=interpolation, cmap=cmap[i])
        # axs.set(title=list(data_dict)[i])
        values = np.unique(data[...,i].ravel()).astype(int)
        colors = [ im.cmap(im.norm(value)) for value in values]
        # create a patch (proxy artist) for every color 
        patches = [ mpatches.Patch(color=colors[j], label=f"${values[j]}$") for j in range(len(values)) ]
        # put those patched as legend-handles into the legend
        axs[i].legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. , ncol=2-i, title=rf'{list(data_dict)[i]}')
        

    # axes labels:
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    fig.set_size_inches(w=LATEX_WIDTH, h=8) 
    plt.savefig(title+'.pdf', bbox_inches = "tight")

    label=[]
    label.append(r'Building the phase diagram for a spin-rtiplet SSH model. Above: seeds are perculated evenly spaced throughout phase space. At each location, two seeds are placed, corresponding to the unitary and non-unitary initial mean-fields. Below: from these seeds (refered to as Layer 0), layers are grown outwards from their parent seed. The outer layers of various seeds overlap, and that with the minimum free energy is selected.')
    text=' '.join(label)
    with open(title+'.txt', 'w') as f:
        f.write(text)

plot_perculation()
