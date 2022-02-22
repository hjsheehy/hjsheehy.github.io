from lib import *
V=10.28
w=1
vv=np.arange(0,4,0.02)
xmin,xmax,ymin,ymax=min(vv),max(vv),-2,2
le,lv=401,len(vv)
mat=np.zeros([lv,le])
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
    mu=0

    n_cells=43

    tb.cut_piece(n_cells, [0], glue_edgs=True)
    tb.set_onsite(-mu,orbital='s')
    tb.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B', label='v')
    tb.set_hopping(-w,hop_vector=[1,0],atom_i='B',atom_f='A', label='w')
    tb.add_impurities(V,[0,0],label='V')
    #########################################################

    tb.solve()

    energy_interval=np.linspace(ymin,ymax,le)
    resolution=0.1
    tb.calculate_greens_function(energy_interval,resolution)

    model=tb
    
    for j,energy in enumerate(energy_interval):
        ldos=tb.local_density_of_states(energy=energy, atom=None, orbital=None)
        dos=np.mean(ldos)

        mat[i,j]=dos

    ##################################################

fig, axs = plt.subplots(1, 1, sharex='all', sharey='all')

im=axs.imshow(mat.T,extent=[xmin,xmax,ymin,ymax], origin='lower')
plt.title(r'Normal SSH')
divider = make_axes_locatable(axs)
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar=fig.colorbar(im,cax=cax)
axs.set_xlabel(r'$v/w$')
axs.set_ylabel(rf'$\overline{{\hat{{\mathcal{{G}}}}(\omega)}}$')

fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
title=r'normal-SSH-phase-diagram'
title=os.path.join(FIG,title)
plt.savefig(title+'.pdf', bbox_inches = "tight")
text=rf'Phase diagram of the SSH model in the normal state with $\mu/w={mu}$. The density of states are plotted as a function of the applied bias energy $\omega$ and intracell, interorbital hopping parameter $v/w$. The topological transition occurs at $v/w=1$. The effect of the chemical potential is to shift the phase boundary, whereas the impurity shifts the phase boundary on a single site, which we observe as a faint line in the density of states, shifted by $V={V}$.'
with open(title+'.txt', 'w') as f:
    f.write(text)

fig, axs = plt.subplots(1, 1, sharex='all', sharey='all')
fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
fig,axs=model.plot_unit_cell(fig,axs,s=150)
title=r'SSH-unit-cell'
title=os.path.join(FIG,title)
plt.savefig(title+'.pdf', bbox_inches = "tight")
text=rf'The dimer unit cell of the SSH model is depicted, together with hopping parameters $w$ and $v$, as well as the lattice basis $b_0$ and $b_1$. Atoms $A$ and $B$ are colored in blue and orange, respectively.'
with open(title+'.txt', 'w') as f:
    f.write(text)
