from lib import *
V=20.28
w=1
fig, axs = plt.subplots(4, 1, sharex='all', sharey='all')
v=0.67
muu=[0,-0.57]
for i,mu in enumerate(muu):
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
    x=11
    tb.add_impurities(V,[-x,0],label='V')
    tb.add_impurities(V,[x,0],label='V')
    #########################################################

    tb.solve()

    energy_interval=np.linspace(-4,+4,401)
    resolution=0.1
    tb.calculate_greens_function(energy_interval,resolution)

    fig, axs[0+2*i] = tb.plot_lattice(fig, axs[0+2*i], energy=0, atoms=None, plot_ldos=True, plot_magnetism=False, s=10)
    fig, axs[1+2*i] = tb.plot_lattice(fig, axs[1+2*i], energy=0, atoms=None, plot_ldos=False, plot_magnetism=True, s=10)
    label=rf'$v/w={v}, \mu/w={muu[i]}$'   
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

title=rf'normal-SSH-LDOS-two-impurities'
title=os.path.join(FIG,title)
plt.savefig(title+'.pdf', bbox_inches = "tight")
text=rf'The local density of states and magnetism of the SSH model in the normal state with two impurites with strong coupling strengths $V={V}$ at $x=\pm11$. Curiously, a kink forms in the middle $x=0$, where there is no impurity.'

with open(title+'.txt', 'w') as f:
    f.write(text)
