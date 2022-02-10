from plt_lib import *
########################################################
########################## Data ########################
########################################################
xx=[]
yy=[]
zz1=[]
zz2=[]
zz3=[]
zz4=[]
zz5=[]
zz6=[]
data_points=[]

plot_ldos=False

for filename in filenames:
    globals().update(conf_file(filename))

    import _pickle as cPickle
    with open(filename, 'rb') as f:
        data = cPickle.load(f)


    ########################################################################################
    ######################## 
    ########################################################################################
    if plot_ldos:
        n_pts=40    

        # fig,axs = plot_data.lattice(energy=energy, model=data, layer=layer, spins=np.arange(plot_data.n_spins), orbitals=np.arange(data.n_orbitals), annotate_orbs=True, show_cell_borders=True, show_basis_vectors=False, show_hopping=True, show_impurities=True, show_only_centre=True)

        fig, axs = plt.subplots(3,figsize=(LATEX_WIDTH, 7.5), gridspec_kw={'width_ratios':[1], 'height_ratios':[20,20,1]})

        energy=0
        dimension=0

        plt_data = plot_data(fig,axs[2],data)

        fig, axs[2] = plt_data.differential_current_map(energy, layer=(ALL,ALL), orbital=[0,1], spin=None, cartesian=True, n_pts=n_pts, gaussian_mean=0.1)

        label_x = r'$x/a_x$'
        label_y = r'$y/a_y$'
        axs[2].set_ylabel(label_y)
        axs[2].set_xlabel(label_x)

        axs[2].set_title('Zero-bias local density of states')

        plt_data = plot_data(fig,axs[1],data)

        fig,axs[1] = plt_data.band_structure(dimension)
        axs[1].set_title('Bandstructure')
        label_x = r''
        axs[1].set_xlabel(label_x)

        colors = [cm.rainbow(i) for i in np.linspace(0, 1, 4)]

        centre=data.centre
        ldos_a = data.local_density_of_states(energy, orbital=0)
        ldos_b = data.local_density_of_states(energy, orbital=1)
        stag=ldos_a-ldos_b
        ldos=ldos_a+ldos_b
        stag=stag[:,centre[1]]
        ldos=ldos[:,centre[1]]
        vmin,vmax=np.amin(ldos),np.amax(ldos)
        gorkov=data.gorkov()

        y=centre[1]
        Delta_T=np.real([[gorkov[x,y,0,s,x,y,1,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
        Delta_T=np.sum(Delta_T,0)/data.n_spins
        Delta_R=np.real([[gorkov[x-1,y,1,s,x,y,0,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
        Delta_R=np.sum(Delta_R,0)/data.n_spins

        x=np.arange(-centre[0],centre[0]+1)
        for k,y in enumerate([ldos,stag,Delta_T,Delta_R]):
            y=np.fft.fftshift(y)
            axs[0].plot(x,y, 
                    color=colors[k])
        axs[0].legend([r'Zero-bias LDOS',r'Staggered density', r'$\Delta_T$',r'$\Delta_R$'],labelcolor='k',loc='upper right')

        label_x = r''
        label_y = r'Field strength'
        axs[0].set_ylabel(label_y)
        axs[0].set_xlabel(label_x)

        axs[0].set_title(f'$U_T={U_T}$, $U_R={U_R}$, $\mu={mu}, V={V}$, {state}')

        plt.tight_layout()

        title=f'U_T={U_T}, U_R={U_R}, mu={mu}, V={V}, {state}, {trajectory}'
        plt.savefig(output+title+'_.pdf', bbox_inches = "tight")

        plt.close()

    ########################################################################################


    xx.append(U_T)
    yy.append(U_R)

    energy=0

    centre=data.centre
    ldos_a = data.local_density_of_states(energy, orbital=0)
    ldos_b = data.local_density_of_states(energy, orbital=1)
    stag_den=ldos_a-ldos_b
    stag_den=stag_den[:,centre[1]]
    l=len(stag_den)
    z1=np.sum(np.abs(stag_den))/l
    z4=np.sum(stag_den)**2/np.sum(np.square(stag_den))
    M0=data.magnetism(energy,orbital=0) 
    M1=data.magnetism(energy,orbital=1)
    M=M0+M1
    z5=np.abs(np.sum(M)/l)
    gorkov=data.gorkov()

    y=centre[1]
    # Delta_T=np.abs([[gorkov[x,y,0,s,x,y,1,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
    Delta_T=np.abs([gorkov[x,y,0,1,x,y,1,1] for x in range(data.dimensions[0])])
    # l=len(Delta_T[0])
    # if state=='INT':
    #     Delta_T=np.sum(Delta_T)/l
    # if state=='unitary':
    #     Delta_T=np.sum(Delta_T)/(2*l)
    Delta_T=np.sum(Delta_T)/l
    z2=Delta_T
    Delta_R=np.abs([[gorkov[x-1,y,1,s,x,y,0,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
    Delta_R=np.sum(Delta_R)/(2*l) #2 spins
    z3=Delta_R

    z6=data.free_energy[-1]

    zz1.append(z1)
    zz2.append(z2)
    zz3.append(z3)
    zz4.append(z4)
    zz5.append(z5)
    zz6.append(z6)
    data_points.append([state,V,mu,trajectory])

del y

titles=[r'Staggered density',r'$\overline{|\Delta_T|}$',r'$\overline{|\Delta_R|}$', r'Inverse participation ratio (IPR)',r'Magnetisation, $|\overline{\langle\hat M \rangle}|$',r'Free energy']
data_points=np.array(data_points)
data_pointss=np.unique(data_points,axis=0)
for data_point in data_pointss:
    [state,V,mu,trajectory]=data_point
    for i,zz in enumerate([zz1,zz2,zz3,zz4,zz5,zz6]):
        x=np.array(xx); y=np.array(yy); z=np.array(zz)
        
        x=x[np.all(data_points==data_point,axis=1)]
        y=y[np.all(data_points==data_point,axis=1)]
        z=z[np.all(data_points==data_point,axis=1)]
        fig, axs = plt.subplots(1)
        title=titles[i]
        title=title+f', $\mu={mu}, V={V}$, {state}, {trajectory}'

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

        # im=plt.scatter(x, y, c=z, s=42, marker='s')
        im=plt.scatter(x, y, c=z, s=160, marker='s')

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

