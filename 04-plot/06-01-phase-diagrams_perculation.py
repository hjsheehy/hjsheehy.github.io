from plt_lib import *

FILENAMES=matrix_filenames()
xlabel=r'$U_v$'
ylabel=r'$U_w$'
interpolation='none'

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

plot_perculation()
