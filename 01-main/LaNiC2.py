from lib import *

filename=sys.argv[0].split('.')[0]
DATA=os.path.join(DATA,filename)
FIG=os.path.join(FIG,filename)
for directory in [FIG]:
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
    
    tb.cut(n_cells, axes=0, glue_edgs=glue_edgs)
    tb.cut(n_cells, axes=1, glue_edgs=True)
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

    # with open(DATA+'.npz', 'wb') as f:
    #     cPickle.dump([greens_function_kq, tb], f)
    return greens_function_kq, tb

# Plotting:
    
def unit_cell(model):
    FIGNAME='unit_cell'

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    plt.tight_layout()
    model.plot_unit_cell(fig, ax, atoms='all', s=100)
    OUTPUT=os.path.join(FIG,FIGNAME)
    plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight")
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
    plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight")

    with open(OUTPUT+'.txt', 'w') as f:
        f.write(rf'''The local density of states of a Weyl SSH model, that is, a spinless square lattice
with ${n_cells}\times{n_cells}$
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
    plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight")

    with open(OUTPUT+'.txt', 'w') as f:
        f.write(rf'''The spectrum of the simple $\text{{NiC}}_2$ model in the normal state, in k-space, integrated over the direction along the SSH chains and summed over the atomic sites. The SSH dimers are in their trivial state with $\mu={mu}, v={v}, w={w}, t_d={td}$. Curiously, Majorana modes remain, as well as the Bogoliubov-Fermi arc, with a $k=\pi$ phase difference.
''')

def k_space_topological(greens_function):
    FIGNAME='k_space_topological'

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, axes=['integrated','resolved'],omega_min=-2,omega_max=2,vmin=0,vmax=80, atom='integrated')

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH/2) 
    plt.tight_layout()
    OUTPUT=os.path.join(FIG,FIGNAME)
    plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight")

    with open(OUTPUT+'.txt', 'w') as f:
        f.write(rf'''The spectrum of the simple $\text{{NiC}}_2$ model in the normal state, in k-space, integrated over the direction along the SSH chains and summed over the atomic sites. The SSH dimers are in their topological state with $\mu={mu}, v={v}, w={w}, t_d={td}$.
Dirac nodes are observed, which pin Majorana zero-modes.
Between the Dirac nodes is a line of zero-energy states, known as a Bogoliubov-Fermi arc.
''')

def k_space_phase_diagram(calculate):
    
    fignames=['phase_diagram_closed','phase_diagram_open']
    glue_edgss=[true,false]

    for i in range(2):
        figname=fignames[i]
        glue_edgs=glue_edgss[i]

        r=c=4
        vv=[0,0.6,1.1,1.2,1.3,2.4]
        tdd=[0,0.6,1.1,1.2,1.3,2.4]
        
        if calculate:
            greens_list=[]
            for td in tdd:
                green_list=[]
                for v in vv:
                    greens_function_kq, tb = main(v,td,glue_edgs)
                    green_list.append(greens_function_kq)
                greens_list.append(green_list)

            with open(data+figname+'.npz', 'wb') as f:
                cpickle.dump([tdd,vv,greens_list], f)
        else:
            [tdd,vv,greens_list]=np.load(data+figname+'.npz', allow_pickle=true)

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
                        axs[i,j].set_ylabel(rf'$v/w={v}$')
                if i==0:
                    axs[i,j].xaxis.set_label_position("top")
                    axs[i,j].set_xlabel(rf'${td}$')
                    if j==0:
                        axs[i,j].set_xlabel(rf'$t_d={td}$')
                        

        fig.suptitle(greens_function.title)
        fig.supxlabel(greens_function.xlabel)
        fig.supylabel(greens_function.ylabel)
        fig.set_size_inches(w=latex_width, h=1.2*latex_width) 
        plt.subplots_adjust(wspace=0.3, hspace=0.25)
        # plt.tight_layout()

        output=os.path.join(fig,figname)
        plt.savefig(output+'.pdf', bbox_inches = "tight")
        
        if glue_edgs:
            with open(output+'.txt', 'w') as f:
                f.write(rf'''the spectrum of the simple $\text{{nic}}_2$ model in the normal state, in k-space, integrated over the direction along the ssh chains and summed over the atomic sites.
the model consists of an ${n_cells}\times{n_cells}$ lattice, with closed boundary conditions.
the chemical potential is $\mu={mu}$, the interdimer hopping is fixed $w={w}$, and the intradimer $v$ and interchain $t_d$ hoppings are varied.'''+r'''
compare with the model with open boundary conditions fig \ref{fig:phase_diagram_open}.
''')
        else:
            with open(output+'.txt', 'w') as f:
                f.write(r'''the model and quantities depicted as in fig \ref{fig:phase_diagram_closed}, but with the boundary conditions closed.
we observe the diagonal hopping $t_d$ has the effect of {\it squeezing} the topological modes into a narrower interval $\mathcal{i}\subset[-\pi,\pi]$. 
at the conventional ssh topological transition point $v=w$, the narrow interval transitions to its complement interval $\mathcal{i}^c$.
notice the universality with the $\text{nic}_2$ model fig \ref{fig:weyl_phase_diagram_open}, \ref{fig:weyl_phase_diagram_closed}.
''')

#############################################################################
################################# Main ######################################
#############################################################################
mu=0.0
s=0.0
# td=0.9
# v=0.6
w=1.2
n_cells=43

# greens_function_kq, tb = main()

# [greens_function_kq, bdg] = np.load(DATA, allow_pickle=True)

CALCULATE=False
k_space_phase_diagram(CALCULATE)
exit()
# unit_cell(tb)
# ldos_each_atom(greens_function_kq)
k_space_trivial(greens_function_kq)

v=1.2
w=0.6
greens_function_kq, tb = main()
k_space_topological(greens_function_kq)
