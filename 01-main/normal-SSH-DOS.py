from lib import *
V=20.28
w=1
vv=np.arange(0,4,0.02)
xmin,xmax,ymin,ymax=min(vv),max(vv),-2,2
le,lv=401,len(vv)
mat=np.zeros([4,lv,le])
for i,v in enumerate(vv):
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
    
    for j,energy in enumerate(energy_interval):
        ldos=tb.local_density_of_states(energy=energy, atom=None, orbital=None)
        dos=np.mean(ldos)
        mat[0,i,j]=dos
        stag_den=tb.mean_abs_staggered_density(atom_i='A', atom_f='B', energy=energy)
        mat[1,i,j]=stag_den
        stag_ipr=tb.IPR_abs_staggered_density(atom_i='A', atom_f='B', energy=energy)
        mat[2,i,j]=stag_ipr
        mag=tb.mean_abs_magnetism(energy=energy, atom=None, orbital=None)
        mat[3,i,j]=mag

    ##################################################

fig, axs = plt.subplots(2, 2, sharex='all', sharey='all')

im=axs[0,0].imshow(mat[0].T,extent=[xmin,xmax,ymin,ymax], origin='lower', cmap='Blues')
axs[0,0].set_title(r'Normal SSH $\overline{{\hat{{\mathcal{{G}}}}(\omega-\mu)}}$')
divider = make_axes_locatable(axs[0,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar=fig.colorbar(im,cax=cax)

im=axs[0,1].imshow(mat[3].T,extent=[xmin,xmax,ymin,ymax], origin='lower',cmap='Reds')
axs[0,1].set_title(r'$\overline{\langle|\hat{M}|\rangle}$')
divider = make_axes_locatable(axs[0,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar=fig.colorbar(im,cax=cax)

im=axs[1,0].imshow(mat[1].T,extent=[xmin,xmax,ymin,ymax], origin='lower', cmap='Greens')
axs[1,0].set_title(r'$\overline{|\text{Staggered density}|}$')
divider = make_axes_locatable(axs[1,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar=fig.colorbar(im,cax=cax)

im=axs[1,1].imshow(mat[2].T,extent=[xmin,xmax,ymin,ymax], origin='lower', cmap='RdPu')
axs[1,1].set_title(r'$\text{IPR}(|\text{Staggered density}|)$')
divider = make_axes_locatable(axs[1,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar=fig.colorbar(im,cax=cax)


fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
xlabel=r'$v/w$'
ylabel=r'$\omega-\mu$'
plt.xlabel(xlabel)
plt.ylabel(ylabel)

fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
title=r'normal-SSH-phase-diagram'
title=os.path.join(FIG,title)
plt.savefig(title+'.pdf', bbox_inches = "tight")
text=rf'Phase diagram of the SSH model with $n={n_cells}$ dimers and periodic boundary conditions. The density of states are plotted as a function of the applied bias energy $\omega-\mu$, noting that the chemical potential contributes only an overall shift, and as a function of the intracell, interorbital hopping parameter $v/w$. The topological transition occurs at $v/w=1$ and notice the faint density at $\omega=0$ below the transition point, due to the edge states. A large impurity coupling $V={V}$ gaps out any modes on the middle site, effectively cutting the closed chain into open boundary conditions. In order to have a measure of the topological modes, we plot the mean absolute magnetism and the mean absolute staggered density, which signal charge density waves. The inverse participation of the staggered density is a measure of long-range order, demonstrating that the topological modes are localised around the impurity.'
with open(title+'.txt', 'w') as f:
    f.write(text)

fig, axs = plt.subplots(1, 1, sharex='all', sharey='all')
fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
fig,axs=tb.plot_unit_cell(fig,axs,s=150)
title=r'SSH-unit-cell'
title=os.path.join(FIG,title)
plt.savefig(title+'.pdf', bbox_inches = "tight")
text=rf'The dimer unit cell of the SSH model is depicted, together with hopping parameters $w$ and $v$, as well as the lattice basis $b_0$ and $b_1$. Atoms $A$ and $B$ are colored in blue and orange, respectively.'
with open(title+'.txt', 'w') as f:
    f.write(text)
