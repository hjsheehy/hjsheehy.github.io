from plt_lib import *
########################################################
########################## Data ########################
########################################################
xx=[]
yy=[]
zz1=[]
zz2=[]
zz3=[]
MU=[]
n_data=len(filenames)

for k in range(n_data):
    filename = filenames[k]
    globals().update(conf_file(filename))

    import _pickle as cPickle
    with open(filename, 'rb') as f:
        data = cPickle.load(f)
    xx.append(U_T)
    yy.append(U_R)

    energy=0

    centre=data.centre
    ldos_a = data.local_density_of_states(energy, orbital=0)
    ldos_b = data.local_density_of_states(energy, orbital=1)
    ldos=ldos_a-ldos_b
    ldos=ldos[:,centre[1]]
    z1=np.sum(np.abs(ldos))/len(ldos)
    gorkov=data.gorkov()

    y=centre[1]
    Delta_T=np.real([[gorkov[x,y,0,s,x,y,1,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
    Delta_T=np.sum(Delta_T,0)/data.n_spins
    z2=np.sum(Delta_T)/len(Delta_T)
    Delta_R=np.real([[gorkov[x-1,y,1,s,x,y,0,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
    Delta_R=np.sum(Delta_R,0)/data.n_spins
    z3=np.sum(Delta_R)/len(Delta_R)

    zz1.append(z1)
    zz2.append(z2)
    zz3.append(z3)
    MU.append(mu)

del y

titles=['CDW','$\Delta_T$','$\Delta_R$']
MUs=np.unique(MU)
for mu in MUs:
    for i,zz in enumerate([zz1,zz2,zz3]):
        x=np.array(xx); y=np.array(yy); z=np.array(zz)
        
        x=x[MU==mu]
        y=y[MU==mu]
        z=z[MU==mu]
        fig, axs = plt.subplots(1)
        title=titles[i]
        title=title+f', $\mu={mu}$'

        # eps=0.0001
        # z[z>eps]=1

        # Set up a regular grid of interpolation points
        # xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
        # xi, yi = np.meshgrid(xi, yi)

        # # Interpolate
        # from scipy import interpolate
        # rbf = interpolate.Rbf(x, y, z, function='linear')
        # zi = rbf(xi, yi)

        # im=axs.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
        #            extent=[x.min(), x.max(), y.min(), y.max()])
        #plt.scatter(x, y, c=z)
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

        im=plt.scatter(x, y, c=z, s=42, marker='s')

        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.88, 0.15, 0.04, 0.6])
        fig.colorbar(im, cax=cbar_ax)

        label_x = r'$U_T$'
        label_y = r'$U_R$'
        axs.set_ylabel(label_y)
        axs.set_xlabel(label_x)
        axs.set_title(title)

        plt.savefig(output+f'_{title}.pdf', bbox_inches = "tight")
        plt.close()
