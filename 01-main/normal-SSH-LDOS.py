from lib import *
V=10.28
w=1
vv=np.arange(0,4,0.02)
# vv=np.arange(0,4,2.1)
xmin,xmax,ymin,ymax=min(vv),max(vv),-4,4
le,lv=401,len(vv)
mat=np.zeros([lv,le])
phase=[r'Trivial',r'Topological']
fig, axs = plt.subplots(4, 1, sharex='all', sharey='all')
vv=[1.39,0.67]
muu=[0,-0.57]
MU=['zero','nonzero']
for j,mu in enumerate(muu):
    for i,v in enumerate(vv):
        #########################################################
        ################# Simple square tb ###################
        #########################################################
        A=Atom([0,0],'A')
        B=Atom([0.25,1],'B')
        A.add_orbital('s')
        B.add_orbital('s')
        lattice_vectors=[[1,0],[0,1]]
        tb=Tightbinding(lattice_vectors,'SSH')
        tb.add_atom(A)
        tb.add_atom(B)
        tb.n_spins=2

        n_cells=43

        tb.cut_piece(n_cells, [0], glue_edgs=True)
        tb.set_onsite(-mu,orbital='s')
        tb.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B', label='v')
        tb.set_hopping(-w,hop_vector=[1,0],atom_i='B',atom_f='A', label='w')
        tb.add_impurities(V,[0,0],label='V')
        #########################################################

        tb.solve()

        energy_interval=np.linspace(xmin,xmax,le)
        resolution=0.1
        tb.calculate_greens_function(energy_interval,resolution)

        model=tb

        fig, axs[0+2*i] = model.plot_lattice(fig, axs[0+2*i], energy=0, atoms=None, plot_ldos=True, plot_magnetism=False, s=10)
        fig, axs[1+2*i] = model.plot_lattice(fig, axs[1+2*i], energy=0, atoms=None, plot_ldos=False, plot_magnetism=True, s=10)
        label=phase[i]+': '+rf'$v/w={vv[i]}$'   
        axs[2*i].annotate(label,
                    xy=(0, 1), xycoords='axes fraction',
                    xytext=(2, -2), textcoords='offset pixels',
                    horizontalalignment='left',
                    verticalalignment='top')

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    xlabel=r'$\hat{x}$'
    ylabel=r'$\hat{y}$'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    title=rf'normal-SSH-LDOS-'+MU[j]
    title=os.path.join(FIG,title)
    plt.savefig(title+'.pdf', bbox_inches = "tight")
    if j==0:
        text=rf'The local density of states and magnetism of the SSH model in the normal state with chemical potential $\mu={mu}$ and strong impurity coupling strength $V={V}$ on the middle site, in order to cut the chain without weak coupling effects. The dark lines represent hopping bonds and the sites are denoted by circles with color corresponding to density. Importantly, we note non-trivial orbital and magnetic characteristics to the edge modes.'
    if j==1:
        text=rf'The local density of states and magnetism of the SSH model in the normal state with chemical potential $\mu={mu}$ and strong impurity coupling strength $V={V}$ on the middle site. Away from zero chemical potential, we observe a charge density wave in the topological phase.'

    with open(title+'.txt', 'w') as f:
        f.write(text)
