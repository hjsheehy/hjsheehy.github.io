from lib import *
V=20.28
w=1
v=1.45


fig, axs = plt.subplots(2, 1, sharex='all', sharey='all')

for i,topology in enumerate(['Trivial','Topological']):

    A=Atom([0,0],'A')
    A.add_orbital('s')
    if topology=='Trivial':
        B=Atom([0.25,0.4],'B')
    if topology=='Topological':
        B=Atom([0.75,0.4],'B')
    B.add_orbital('s')
    lattice_vectors=[[1,0],[0,0.5]]
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
    fig,axs[i]=tb.plot_unit_cell(fig,axs[i],s=150)
    txt=axs[i].annotate(text=rf'{topology} phase', xytext=[0.05,0.43],xy=[0.05,0.35])
    txt.set_bbox(dict(facecolor='w', alpha=0.7, edgecolor='k') )

    axs[i].set_xlabel(r'')
    axs[i].set_ylabel(r'')

fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
xlabel=r'$\hat{x}$'
ylabel=r'$\hat{y}$'
plt.xlabel(xlabel)
plt.ylabel(ylabel)
fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
title=r'SSH-unit-cell'
title=os.path.join(FIG,title)
plt.savefig(title+'.pdf', bbox_inches = "tight")
text=rf'The dimer unit cell of the SSH model depicted together with hopping parameters $w$ and $v$, lattice basis $b_0$ and $b_1$ and atoms $A$ (blue) and $B$ (orange)). Above: trivial dimer phase, in which $v>w$. Below: $v<w$, model exhibits topologically protected edge modes. Our unconventional model includes the possibility of intracell multiorbital, equal-spin attraction $U_v$ and intercell (Coulomb) repulsion $U_w$.'
with open(title+'.txt', 'w') as f:
    f.write(text)
