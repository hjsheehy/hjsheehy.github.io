from lib import *

FILENAME=sys.argv[0].split('.')[0]
DATA=os.path.join(DATA,FILENAME)
FIG=os.path.join(FIG,FILENAME)
for directory in [FIG,DATA]:
    if not os.path.exists(directory):
        os.makedirs(directory)

def main():
    A=Atom([0,0],'A')
    B=Atom([0.5,0],'B')
    A.add_orbital('s')
    B.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    bdg=BogoliubovdeGennes(lattice_vectors,'2D-Weyl-SSH')
    bdg.add_atom(A)
    bdg.add_atom(B)
    bdg.n_spins=1
    bdg.cut(nx, axes=0, glue_edgs=False)
    bdg.cut(ny, axes=1, glue_edgs=True)
    bdg.set_onsite(-mu+s,atom='A')
    bdg.set_onsite(-mu-s,atom='B')

    bdg.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B',label='$v$')
    bdg.set_hopping(-w,hop_vector=[-1,0],atom_i='A',atom_f='B',label='$w$')
    bdg.set_hopping(-td,hop_vector=[0,1],atom_i='B',atom_f='A',label='$t_d$')
    bdg.set_hopping(-td,hop_vector=[0,-1],atom_i='B',atom_f='A',label='$t_d$')
    bdg.set_hopping(-td,hop_vector=[1,1],atom_i='B',atom_f='A',label='$t_d$')
    bdg.set_hopping(-td,hop_vector=[1,-1],atom_i='B',atom_f='A',label='$t_d$')

    # impurity_wall = [[0,i] for i in range(n_cells)]
    # bdg.add_impurities(V,impurity_wall)

    bdg.add_impurities(V,[0,0])

    # bdg.set_hartree(rho)
    # bdg.set_fock(phi,atom_i='A',atom_f='B')
    # bdg.set_fock(phi,atom_i='B',atom_f='A',hop_vector=[1,0])
    # bdg.set_gorkov(chi,atom_i='A',atom_f='B')
    # bdg.set_gorkov(chi,atom_i='B',atom_f='A',hop_vector=[1,0])

    bdg.set_hubbard_u(-Uv,atom_i='A',atom_f='B',hop_vector=[0,0],orbital_i='s',orbital_f='s')
    # bdg.set_hubbard_u(Uw,atom_i='B',atom_f='A',hop_vector=[1,0])

    # bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B', orbital_i='s', orbital_f='s', spin_i=0, spin_f=0,_print=False)
    # bdg.record_gorkov(location_i=[0,0], location_f=[1,0], atom_i='B', atom_f='A', orbital_i='s', orbital_f='s', spin_i=0, spin_f=0,_print=False)

    # bdg.self_consistent_calculation(friction=0., max_iterations=2, absolute_convergence_factor=0.00001)
    bdg.solve()

    energy_interval=np.linspace(-4,4,601)
    resolution=0.05

    greens_function_xy=GreensFunction(bdg,energy_interval,resolution, k_axes=None)

    greens_function_xq=GreensFunction(bdg,energy_interval,resolution, k_axes=[1])

    greens_function_kq=GreensFunction(bdg,energy_interval,resolution, k_axes=[0,1])

    del bdg.eigenvectors
    del bdg.eigenvalues

    return greens_function_xy, greens_function_xq, greens_function_kq, bdg

def main_tb(v,td,glue_edgs):
    A=Atom([0,0],'A')
    B=Atom([0.5,0],'B')
    A.add_orbital('s')
    B.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    tb=TightBinding(lattice_vectors,'2D-Weyl-SSH')
    tb.add_atom(A)
    tb.add_atom(B)
    tb.n_spins=1
    tb.cut(n_cells, axes=0, glue_edgs=glue_edgs)
    tb.cut(n_cells, axes=1, glue_edgs=True)
    tb.set_onsite(-mu+s,atom='A')
    tb.set_onsite(-mu-s,atom='B')

    tb.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B',label='$v$')
    tb.set_hopping(-w,hop_vector=[-1,0],atom_i='A',atom_f='B',label='$w$')
    tb.set_hopping(-td,hop_vector=[0,1],atom_i='B',atom_f='A',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[0,-1],atom_i='B',atom_f='A',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[1,1],atom_i='B',atom_f='A',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[1,-1],atom_i='B',atom_f='A',label='$t_d$')

    tb.solve()

    energy_interval=np.linspace(-4,4,601)
    resolution=0.05

    greens_function_kq=GreensFunction(tb,energy_interval,resolution, k_axes=[0,1])

    del tb.eigenvectors
    del tb.eigenvalues

    return greens_function_kq

# Plotting:
    
def unit_cell(model):
    FIGNAME='unit_cell'

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    plt.tight_layout()
    model.plot_unit_cell(fig, ax, atoms='all', s=100)
    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight", dpi=DPI)
    # plt.close()

    with open(output+'.txt', 'w') as f:
        f.write(rf'''Two-dimensional multiorbital Weyl SSH model.
Along the $\hat{{x}}$-direction, atomic sites A and B within the unit cell form an SSH chain, with intracell hopping $v$ and intercell hopping $w$. 
An additional interchain hopping parameter couples the chains into a two-dimensional stack. 
Following Rosenberg and Manousakis, we study a diagonal hopping $t_d$, which is a Type II extension of the SSH model in the classification according to Bo-Hung [CITE: Two-Dimensional Extended
Su-Schrieffer-Heeger Model; Masters Thesis; Chen, Bo-Hung; July 2018].
The basis for the lattice is given by the vectors $b_0$ and $b_1$.''')

def ldos_each_atom(greens_function):
    FIGNAME='ldos_each_atom'
    omega=0

    fig, ax = plt.subplots(1, 2, sharey='row')
    
    bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
    n=int(n_cells/2-3)
    for i,atom in enumerate(['A','B']):
        ax[i]=greens_function.plot_ldos(ax[i],energy=0,atom=atom)
        ax[i].set_title('')
        ax[i].set_xlabel('')
        ax[i].set_ylabel('')
        ax[i].text(n, n, "Atom "+atom, ha="right", va="top", size=10,
        bbox=bbox_props)

    fig.suptitle(greens_function.title)
    fig.supxlabel(greens_function.xlabel)
    fig.supylabel(greens_function.ylabel)
    fig.set_size_inches(w=LATEX_WIDTH, h=0.6*LATEX_WIDTH) 
    plt.tight_layout()
    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight", dpi=DPI)

    with open(output+'.txt', 'w') as f:
        f.write(rf'''The local density of states of a Weyl SSH model, that is, a spinless square lattice
with ${nx}\times{ny}$
lattice cells, with chemical potential $\mu/t={mu:.2f}$, zero shift $s=0$ and
applied bias $\omega={omega}$ and 
spectral resolution $\epsilon={greens_function.resolution}$.
The SSH chains are in their topological phase $w={w}>v={v}$, with diagonal hopping
$t_d={td}$, 
which does gives rise to no topological modes perpendicular to the chains.
Notice the zero bias density at the left edge for atom A and on the right for atom B, a signature of topological modes.
''')

def real_space(greens_function):
    FIGNAME='ldos'

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, axes=['resolved','integrated'],omega_min=-2,omega_max=2,vmin=0,vmax=80)

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH/2) 
    plt.tight_layout()
    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight", dpi=DPI)

    with open(output+'.txt', 'w') as f:
        f.write(rf'''The spectrum of the Weyl-SSH model, in its normal state, integrated over the direction perpendicular to the SSH chains and summed over the atomic sites within the lattice cell.
Zero-energy modes are observed at the edges, with a faint tail of density.
''')

def k_space(greens_function):
    FIGNAME='k_space'

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, axes=['integrated','resolved'],omega_min=-2,omega_max=2,vmin=0,vmax=80)

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH/2) 
    plt.tight_layout()
    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight", dpi=DPI)

    with open(output+'.txt', 'w') as f:
        f.write(rf'''The spectrum of the Weyl-SSH model in the normal state, in k-space, integrated over the direction along the SSH chains and summed over the atomic sites.
Dirac nodes are observed, which pin Majorana zero-modes.
Between the Dirac nodes is a line of zero-energy states, known as a Bogoliubov-Fermi arc.
''')

def majorana_fermi_arc(greens_function):
    FIGNAME='majorana_fermi_arc'

    fig, ax = plt.subplots()
    
    majorana_energy=0
    Bogoliubov_Fermi_arc_energy=0
    majorana_ky=2.06
    Bogoliubov_Fermi_arc_ky=0.5
    ax = greens_function.plot_spectrum(ax, energy=majorana_energy, axes=['resolved',majorana_ky], omega_min=0,omega_max='default',vmin=0,vmax=6,label='Majorana')
    ax = greens_function.plot_spectrum(ax, energy=Bogoliubov_Fermi_arc_energy, axes=['resolved',Bogoliubov_Fermi_arc_ky], omega_min=0,omega_max='default',vmin=0,vmax=6,label='Bogoliubov-Fermi arc')
    ax.set_title('')
    ax.legend()

    ins = ax.inset_axes([0.2,0.3,0.6,0.4])

    ins = greens_function.plot_spectrum(ins, axes=['integrated','resolved'], omega_min=-1,omega_max=1,vmin='default',vmax='default')
    
    ins.scatter(majorana_ky,majorana_energy,c='r',s=10)
    ins.scatter(Bogoliubov_Fermi_arc_ky,Bogoliubov_Fermi_arc_energy,c='blue',s=20)
    ins.set_xlim([0,np.pi])
    ax.set_title('Quasiparticle spectrum')

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH*0.8) 
    plt.tight_layout()
    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight", dpi=DPI)

    with open(output+'.txt', 'w') as f:
        f.write(rf'''The spectrum of the Majorana and Bogoliubov-Fermi arc quasiparticle modes, at zero-applied bias, with the vertical momentum specified by the (colour-coordinated) spectral plot.
The quasiparticle density is plotted over the direction parallel to the SSH chains. 
The Majorana mode (red) is localised at the edges. 
The Bogoliubov-Fermi arc, in contrast, tails.
''')

def k_space_phase_diagram(CALCULATE):
    
    FIGNAMES=['phase_diagram_closed','phase_diagram_open']
    glue_edgss=[True,False]

    for i in range(2):
        FIGNAME=FIGNAMES[i]
        glue_edgs=glue_edgss[i]

        r=c=4
        vv=[0,0.3,0.6,1.1,1.2,1.3,2.4,5,10,50]
        tdd=[0,0.3,0.6,1.1,1.2,1.3,2.4,5,10,50]
        
        if CALCULATE:
            greens_list=[]
            for td in tdd:
                green_list=[]
                for v in vv:
                    greens_function_kq = main_tb(v,td,glue_edgs)
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
                        axs[i,j].set_ylabel(rf'$v/w={v}$')
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
                f.write(rf'''The spectrum of the simple Weyl-SSH model in the normal state, in k-space, integrated over the direction along the SSH chains and summed over the atomic sites.
The model consists of an ${nx}\times{ny}$ lattice, with closed boundary conditions.
The chemical potential is $\mu={mu}$, the interdimer hopping is fixed $w={w}$, and the intradimer $v$ and interchain $t_d$ hoppings are varied.'''+r'''
Compare with the model with open boundary conditions Fig \ref{fig:phase_diagram_open}, and notice the universality when comparing to the $\text{NiC}_2$ model Fig \ref{fig:phase_diagram_closed}, \ref{fig:phase_diagram_open}.
''')
        else:
            with open(OUTPUT+'.txt', 'w') as f:
                f.write(r'''The model and quantities depicted as in Fig \ref{fig:weyl_phase_diagram_closed}, but with the boundary conditions open.
We observe the diagonal hopping $t_d$ gives rise to two Dirac cones.
In the SSH topological phase, a Bogliubov-Fermi arc connects them through the centre of the Brillouin zone.
Then, in the SSH trivial phase, we find the Fermi arc instead passes through the edge of the Brillouin zone.
''')

def spectrum(greens_function):
    FIGNAME='spectrum'

    fig, ax = plt.subplots()

    ax = greens_function.plot_energy_spectrum(ax, sites='integrated',atom='integrated',xmin=-2,xmax=2,ymin=0,ymax=2000)

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH/2) 
    plt.tight_layout()
    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight", dpi=DPI)

    with open(output+'.txt', 'w') as f:
        f.write(rf'''The spectrum of the Weyl-SSH model in the normal state with open-boundary conditions, summed over the atomic sites. The parameters are $\mu={mu}, v={v}, w={w}$ and $t_d={td}$. We observe the topological mode as a density peak within the band gap.
''')
#############################################################################
################################# Main ######################################
#############################################################################
mu=0.0
s=0.0
td=0.9
v=0.6
w=1.2
Uv=0
# Uw=1
rho=0
phi=0
chi=0
V=0
n_cells=nx=ny=43

CALCULATE=False
if CALCULATE:
    greens_function_xy, greens_function_xq, greens_function_kq, bdg = main()
    with open(os.path.join(DATA,'main.npz'), 'wb') as f:
        cPickle.dump([greens_function_xy, greens_function_xq, greens_function_kq, bdg], f)
else:
    [greens_function_xy, greens_function_xq, greens_function_kq, bdg] = np.load(os.path.join(DATA,'main.npz'), allow_pickle=True)

unit_cell(bdg)
ldos_each_atom(greens_function_xy)
real_space(greens_function_xy)
k_space(greens_function_kq)
majorana_fermi_arc(greens_function_xq)
spectrum(greens_function_kq)
k_space_phase_diagram(CALCULATE=False)
