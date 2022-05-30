from lib import *

FILENAME=sys.argv[0].split('.')[0]
DATA=os.path.join(DATA,FILENAME)
FIG=os.path.join(FIG,FILENAME)
for directory in [FIG,DATA]:
    if not os.path.exists(directory):
        os.makedirs(directory)

def main():
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

    impurity_wall = [[0,i] for i in range(n_cells)]
    bdg.add_impurities(V,impurity_wall)

    # bdg.add_impurities(V,[0,0])

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

# Plotting:
def ldos_each_atom(greens_function):
    FIGNAME='ldos_each_atom'

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
    if caption:
        with open(output+'.txt', 'w') as f:
            f.write(rf'''Local density of states with an impurity at the centre. 
The chemical potential is $\mu/t={mu:.2f}$ and the lattice is ${nx}\times{ny}$
with periodic boundary conditions along the vertical, and open boundary conditions
along the horizontal.
The impurity has coupling strength $V/t={V:.2f}$. ''')

def real_space(greens_function):
    FIGNAME='ldos'

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, axes=['resolved','integrated'],omega_min=-2,omega_max=2,vmin=0,vmax=80)

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH/2) 
    plt.tight_layout()
    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight", dpi=DPI)
    
    if caption:
        with open(output+'.txt', 'w') as f:
            f.write(rf'''Local density of states with an impurity at the centre. 
The chemical potential is $\mu/t={mu:.2f}$ and the lattice is ${nx}\times{ny}$
with periodic boundary conditions perpendicular to the axis of the dimers, and open boundary conditions
along the axis parallel.
The impurity wall runs perpendicular to the axis and through the middle, with coupling strength $V/t={V:.2f}$. 
Additional modes are seen although at such weak coupling, it is not clear whether they are gapped.''')

def k_space(greens_function):
    FIGNAME='k_space'

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, axes=['integrated','resolved'],omega_min=-2,omega_max=2,vmin=0,vmax=40)

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH/2) 
    plt.tight_layout()
    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight", dpi=DPI)
    
    if caption:
        with open(output+'.txt', 'w') as f:
            f.write(rf'''Local density of states integrated over the axis of the SSH dimers and summed over the intracell atomic sites show clearly that there exists an additional gapped mode due to the impurity wall. 
The chemical potential is $\mu/t={mu:.2f}$ and the lattice is ${nx}\times{ny}$
with periodic boundary conditions along the axis perpendicular to the dimers, and open boundary conditions parallel.
The impurity wall has coupling strength $V/t={V:.2f}$.''')

def majorana_fermi_arc(greens_function):
    FIGNAME='majorana_fermi_arc'

    fig, ax = plt.subplots()
    
    majorana_energy=0
    majorana_ky=2.06
    Bogoliubov_Fermi_arc_energy=0
    Bogoliubov_Fermi_arc_ky=0.5
    additional_mode_energy=0
    additional_mode_ky=np.pi
    ax = greens_function.plot_spectrum(ax, energy=majorana_energy, axes=['resolved',majorana_ky], omega_min=0,omega_max='default',vmin=0,vmax=6,label='Majorana')
    ax = greens_function.plot_spectrum(ax, energy=Bogoliubov_Fermi_arc_energy, axes=['resolved',Bogoliubov_Fermi_arc_ky], omega_min=0,omega_max='default',vmin=0,vmax=6,label='Bogoliubov-Fermi arc')
    ax = greens_function.plot_spectrum(ax, energy=additional_mode_energy, axes=['resolved',additional_mode_ky], omega_min=0,omega_max='default',vmin=0,vmax=6,label='Impurity mode')
    ax.set_title('')
    ax.legend()

    ins = ax.inset_axes([0.2,0.3,0.6,0.4])

    ins = greens_function.plot_spectrum(ins, axes=['integrated','resolved'], omega_min=-1,omega_max=1,vmin='default',vmax=50)
    
    ins.scatter(majorana_ky,majorana_energy,c='r',s=10)
    ins.scatter(Bogoliubov_Fermi_arc_ky,Bogoliubov_Fermi_arc_energy,c='blue',s=20)
    ins.scatter(additional_mode_ky,additional_mode_energy,c='g',s=20)
    ins.set_xlim([0,np.pi])
    ax.set_title('Quasiparticle spectrum')

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH*0.8) 
    plt.tight_layout()
    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight", dpi=DPI)
    
    if caption:
        with open(output+'.txt', 'w') as f:
            f.write(rf'''The majorana quasiparticle exhibits some pinning around the impurity wall and an additional arc mode exists at $k_y=\pi$, whose nature is to be investigated.''')
#############################################################################
################################# Main ######################################
#############################################################################
A=Atom([0,0],'A')
B=Atom([0.5,0],'B')
A.add_orbital('s')
B.add_orbital('s')
lattice_vectors=[[1,0],[0,1]]
bdg=BogoliubovdeGennes(lattice_vectors,'2D-Weyl-SSH')
bdg.add_atom(A)
bdg.add_atom(B)
bdg.n_spins=1
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
V=1.21
n_cells=nx=ny=43

# greens_function_xy, greens_function_xq, greens_function_kq, bdg = main()
# with open(os.path.join(DATA,'main.npz'), 'wb') as f:
#     cPickle.dump([greens_function_xy, greens_function_xq, greens_function_kq, bdg], f)

[greens_function_xy, greens_function_xq, greens_function_kq, bdg] = np.load(os.path.join(DATA,'main.npz'), allow_pickle=True)

caption=True
ldos_each_atom(greens_function_xy)
real_space(greens_function_xy)
k_space(greens_function_kq)
majorana_fermi_arc(greens_function_xq)
