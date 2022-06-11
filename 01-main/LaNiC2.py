from lib import *

FILENAME=sys.argv[0].split('.')[0]
DATA=os.path.join(DATA,FILENAME)
FIG=os.path.join(FIG,FILENAME)
for directory in [FIG,DATA]:
    if not os.path.exists(directory):
        os.makedirs(directory)


def main(v,td,glue_edgs):
    A=Atom([0,0],'A')
    B=Atom([0.5,0],'B')
    A.add_orbital('s')
    B.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    tb=TightBinding(lattice_vectors,'LaNiC2')
    tb.add_atom(A)
    tb.add_atom(B)
    tb.n_spins=1
    
    tb.cut(nx, axes=0, glue_edgs=glue_edgs)
    tb.cut(ny, axes=1, glue_edgs=True)
    # tb.set_kpts([n_cells,n_cells])
    tb.set_onsite(-mu+s,atom='A')
    tb.set_onsite(-mu-s,atom='B')

    # tb.add_impurities(2.1,[0,0])

    tb.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B',label='$v$')
    tb.set_hopping(-w,hop_vector=[-1,0],atom_i='A',atom_f='B',label='$w$')
    tb.set_hopping(-td,hop_vector=[0,1],atom_i='B',atom_f='A',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[0,-1],atom_i='B',atom_f='A',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[-1,1],atom_i='A',atom_f='B',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[-1,-1],atom_i='A',atom_f='B',label='$t_d$')

    tb.solve()

    energy_interval=np.linspace(-8,8,601)
    resolution=0.05

    greens_function_kq=GreensFunction(tb,energy_interval,resolution, k_axes=[0,1])

    del tb.eigenvectors
    del tb.eigenvalues

    return greens_function_kq, tb

# Plotting:
    
def unit_cell(model):
    FIGNAME='unit_cell'

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    plt.tight_layout()
    model.plot_unit_cell(fig, ax, atoms='all', s=100)
    OUTPUT=os.path.join(FIG,FIGNAME)
    plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight", dpi=DPI)
    # plt.close()

    with open(OUTPUT+'.txt', 'w') as f:
        f.write(r'''Planar model of the $\text{NiC}_2$ planes of $\text{LaNiC}_2$.
Along the $\hat{{x}}$-direction, atomic sites A and B within the unit cell form an SSH chain, with intracell hopping $v$ and intercell hopping $w$. 
An additional interchain hopping parameter $t_d$ couples the chains into a two-dimensional stack. 
''')

def ldos_each_atom(greens_function):
    FIGNAME='ldos_each_atom'

    fig, ax = plt.subplots(1, 2, sharey='row')

    bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
    n=0.1#int(n_cells/2-3)

    for i,atom in enumerate(['A','B']):
        ldos=greens_function.local_density_of_states(energy='resolved',atom=atom)
        ldos=np.fft.fftshift(ldos)
        ldos=np.argmax(np.abs(ldos),axis=2)
        g=greens_function.energy_interval
        emin,emax,ne=np.min(g),np.max(g),len(g)
        ldos=ldos/ne-emin
        pi=np.pi
        ax[i].imshow(ldos.T,origin='lower',extent=[-pi,pi,-pi,pi])
        ax[i].set_title('')
        ax[i].set_xlabel('')
        ax[i].set_ylabel('')
        ax[i].text(n, n, "atom "+atom, ha="right", va="top", size=10,
        bbox=bbox_props)

    # fig.suptitle(greens_function.title)
    # fig.supxlabel(greens_function.xlabel)
    # fig.supylabel(greens_function.ylabel)
    fig.set_size_inches(w=LATEX_WIDTH, h=0.6*LATEX_WIDTH) 
    plt.tight_layout()
    OUTPUT=os.path.join(FIG,FIGNAME)
    plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight", dpi=DPI)

    with open(OUTPUT+'.txt', 'w') as f:
        f.write(rf'''The local density of states of a Weyl SSH model, that is, a spinless square lattice
with ${nx}\times{ny}$
lattice cells, with chemical potential $\mu/t={mu:.2f}$, zero shift $s=0$ and
spectral resolution $\epsilon={greens_function.resolution}$.
The SSH chains are in their topological phase $w={w}>v={v}$, with diagonal hopping
$t_d={td}$, 
which does gives rise to no topological modes perpendicular to the chains.
Notice the zero bias density at the left edge for atom A and on the right for atom B, a signature of topological modes.
''')

def k_space_trivial(greens_function):
    FIGNAME='k_space_trivial'

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, axes=['integrated','resolved'],omega_min=-2,omega_max=2,vmin=0,vmax=80, atom='integrated')

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH/2) 
    plt.tight_layout()
    OUTPUT=os.path.join(FIG,FIGNAME)
    plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight", dpi=DPI)

    with open(OUTPUT+'.txt', 'w') as f:
        f.write(rf'''The spectrum of the simple $\text{{NiC}}_2$ model in the normal state, in k-space, integrated over the direction along the SSH chains and summed over the atomic sites. 
        The boundary conditions are open along the axis of the dimers and closed in the perpendicular direction.
The SSH dimers are in their trivial state with $\mu={mu}, v={v}, w={w}, t_d={td}$. 
The Bogoliubov-Fermi arc links the Weyl nodes through the outside of the Brillouin zone.
''')

def k_space_topological(greens_function):
    FIGNAME='k_space_topological'

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, axes=['integrated','resolved'],omega_min=-2,omega_max=2,vmin=0,vmax=80, atom='integrated')

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH/2) 
    plt.tight_layout()
    OUTPUT=os.path.join(FIG,FIGNAME)
    plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight", dpi=DPI)

    with open(OUTPUT+'.txt', 'w') as f:
        f.write(rf'''The spectrum of the simple $\text{{NiC}}_2$ model in the normal state, in k-space, integrated over the direction along the SSH chains and summed over the atomic sites. 
The boundary conditions are open along the axis of the dimers and closed in the perpendicular direction.
The SSH dimers are in their topological state with $\mu={mu}, v={v}, w={w}, t_d={td}$.
The Weyl nodes pin Majorana zero-modes, between which runs
a Bogoliubov-Fermi arc through the centre of the Brillouin zone.
''')

def k_space_transition(greens_function):
    FIGNAME='k_space_transition'

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, axes=['integrated','resolved'],omega_min=-2,omega_max=2,vmin=0,vmax=80, atom='integrated')

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH/2) 
    plt.tight_layout()
    OUTPUT=os.path.join(FIG,FIGNAME)
    plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight", dpi=DPI)

    with open(OUTPUT+'.txt', 'w') as f:
        f.write(rf'''The spectrum of the simple $\text{{NiC}}_2$ model in the normal state, in k-space, integrated over the direction along the SSH chains and summed over the atomic sites. 
The boundary conditions are open along the axis of the dimers and closed in the perpendicular direction.
The SSH dimers are at their topological transition point with $\mu={mu}, v=w={w}, t_d={td}$.
''')

def k_space_phase_diagram(CALCULATE):
    
    FIGNAMES=['phase_diagram_closed','phase_diagram_open']
    glue_edgss=[True,False]

    for i in range(2):
        FIGNAME=FIGNAMES[i]
        glue_edgs=glue_edgss[i]

        r=c=4
        vv=[0,0.6,1.1,1.2,1.3,2.4,10,100]
        tdd=[0,0.6,1.1,1.2,1.3,2.4,10,100]
        
        if CALCULATE:
            greens_list=[]
            for td in tdd:
                green_list=[]
                for v in vv:
                    greens_function_kq, tb = main(v,td,glue_edgs)
                    green_list.append(greens_function_kq)
                greens_list.append(green_list)

            with open(os.path.join(DATA,FIGNAME+'.npz'), 'wb') as f:
                cPickle.dump([tdd,vv,greens_list], f)
        else:
            [tdd,vv,greens_list]=np.load(os.path.join(DATA,FIGNAME+'.npz'), allow_pickle=True)

        fig, axs = plt.subplots(len(vv),len(tdd))

        bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
        for i in range(len(vv)):
            for j in range(len(tdd)):
                k=len(vv)-i-1
                greens_function = greens_list[j][k]
                td,v=tdd[j],vv[k]
                axs[i,j] = greens_function.plot_spectrum(axs[i,j], axes=['integrated','resolved'],omega_min=-2,omega_max=2,vmin=0,vmax=20, atom='integrated')
                axs[i,j].set_title('')
                axs[i,j].set_xlabel('')
                axs[i,j].set_ylabel('')
                if i!=len(vv)-1:
                    axs[i,j].set_xticks([])
                if j!=0:
                    axs[i,j].set_yticks([])
                if j==len(tdd)-1:
                    axs[i,j].yaxis.set_label_position("right")
                    axs[i,j].set_ylabel(rf'${v}$')
                    if i==len(vv)-1:
                        axs[i,j].set_ylabel(rf'$v={v}$')
                if i==0:
                    axs[i,j].xaxis.set_label_position("top")
                    axs[i,j].set_xlabel(rf'${td}$')
                    if j==0:
                        axs[i,j].set_xlabel(rf'$t_d={td}$')
                        

        fig.suptitle(greens_function.title)
        fig.supxlabel(greens_function.xlabel)
        fig.supylabel(greens_function.ylabel)
        fig.set_size_inches(w=LATEX_WIDTH, h=1.2*LATEX_WIDTH) 
        plt.subplots_adjust(wspace=0.3, hspace=0.25)
        # plt.tight_layout()

        OUTPUT=os.path.join(FIG,FIGNAME)
        plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight", dpi=DPI)
        
        if glue_edgs:
            with open(OUTPUT+'.txt', 'w') as f:
                f.write(rf'''The spectrum of the simple $\text{{NiC}}_2$ model in the normal state, in k-space, integrated over the direction along the SSH chains and summed over the atomic sites.
The model consists of an ${nx}\times{ny}$ lattice, with closed boundary conditions.
The chemical potential is $\mu={mu}$, the interdimer hopping is fixed $w={w}$, and the intradimer $v$ and interchain $t_d$ hoppings are varied.'''+r'''
Compare with the model with open boundary conditions Figure \ref{fig:phase_diagram_open}.
''')
        else:
            with open(OUTPUT+'.txt', 'w') as f:
                f.write(r'''The model and quantities depicted as in Figure \ref{fig:phase_diagram_closed}, but with the boundary conditions closed.
We observe the diagonal hopping $t_d$ has the effect of creating a boundary in the Brillouin zone, between its centre and its outside, narrowing the centre with increasing $t_d$.
In the conventional SSH model, the topological transition point, $v=w$, is a metallic phase, between a trivial insulator and a Majorana topological insulator.
In this 2D model, the transition point $w=v$ is a topological insulator with a pair of Majorana Dirac cones {\color{red}check this}, which move towards the centre of the Brillouin zone with increasing $t_d$, a second order transition, with the Dirac nodes merging at $t_d\to\infty$.
The inside of the Brillouin zone, between the Dirac cones, exhibits a Fermi arc in the SSH topological phase, $v<w$.
The transition, $v=w$, is a first-order transition; the Fermi arc splits into two, one sharply moving up the spectrum, into the bulk, as the other moves downwards, which is a Majorana topological insulting phase.
The SSH trivial phase, $v>w$, becomes non-trivial in 2D; with the Fermi arc in the opposite side of the Brillouin zone.
Notice the universality with the $\text{NiC}_2$ model Figure \ref{fig:weyl_phase_diagram_open}, \ref{fig:weyl_phase_diagram_closed}.
''')

#############################################################################
################################# Main ######################################
#############################################################################
mu=0.0
s=0.0
td=0.9
v=0.6
w=1.2
n_cells=nx=ny=43

# k_space_phase_diagram(CALCULATE=True)

# glue_edgs=False
# greens_function_kq, tb = main(v,td,glue_edgs)
# with open(os.path.join(DATA,'topological.npz'), 'wb') as f:
#     cPickle.dump([greens_function_kq, tb], f)

[greens_function_kq, tb] = np.load(os.path.join(DATA,'topological.npz'), allow_pickle=True)

unit_cell(tb)
# ldos_each_atom(greens_function_kq)
# k_space_topological(greens_function_kq)

# v=1.2
# w=0.6
# greens_function_kq, tb = main(v,td,glue_edgs)
# with open(os.path.join(DATA,'trivial.npz'), 'wb') as f:
#     cPickle.dump([greens_function_kq, tb], f)

# [greens_function_kq, tb] = np.load(os.path.join(DATA,'trivial.npz'), allow_pickle=True)
# k_space_trivial(greens_function_kq)

# v=1.2
# w=1.2
# greens_function_kq, tb = main(v,td,glue_edgs)
# with open(os.path.join(DATA,'transition.npz'), 'wb') as f:
#     cPickle.dump([greens_function_kq, tb], f)

# [greens_function_kq, tb] = np.load(os.path.join(DATA,'transition.npz'), allow_pickle=True)
# k_space_transition(greens_function_kq)
