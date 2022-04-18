from lib import *

filename=sys.argv[0].split('.')[0]
DATA=os.path.join(DATA,filename+'.npz')
FIG=os.path.join(FIG,filename)
if not os.path.exists(FIG):
    os.makedirs(FIG)

def main():
    bdg.cut(n_cells, axes=0, glue_edgs=False)
    bdg.cut(n_cells, axes=1, glue_edgs=True)
    bdg.set_onsite(-mu+s,atom='A')
    bdg.set_onsite(-mu-s,atom='B')

    bdg.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B',label='$v$')
    bdg.set_hopping(-w,hop_vector=[-1,0],atom_i='A',atom_f='B',label='$w$')
    bdg.set_hopping(-td,hop_vector=[0,1],atom_i='B',atom_f='A',label='$t_d$')
    bdg.set_hopping(-td,hop_vector=[0,-1],atom_i='B',atom_f='A',label='$t_d$')
    bdg.set_hopping(-td,hop_vector=[1,1],atom_i='B',atom_f='A',label='$t_d$')
    bdg.set_hopping(-td,hop_vector=[1,-1],atom_i='B',atom_f='A',label='$t_d$')

    impurity_wall = [[0,i] for i in range(int(n_cells/2))]
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

    with open(DATA, 'wb') as f:
        cPickle.dump([greens_function_xy, greens_function_xq, greens_function_kq, bdg], f)
    return greens_function_xy, greens_function_xq, greens_function_kq, bdg

# Plotting:

def iterations(bdg):
    gorkov_v=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])
    # gorkov_w=bdg.gorkov(atom_i='B', atom_f='A', hop_vector=[1,0])

    plt.plot(np.real(bdg._gorkov_iterations[0]))
    plt.plot(np.real(bdg._gorkov_iterations[1]))
    plt.show()
    plt.close()
    
def unit_cell(model):
    FIGNAME='unit_cell'

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    plt.tight_layout()
    model.plot_unit_cell(fig, ax, atoms='all', s=100)
    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight")
    # plt.close()

    with open(output+'.txt', 'w') as f:
        f.write(rf'''The local density of states of a spinless square lattice tight-binding
model at $\mu/t={mu:.2f}$ with ${nx}\times{ny}$
sites, a single orbital with an impurity at the centre with coupling
strength $V/t={V:.2f}$. The impurity gives rise to Fridel's eponymous
waves in the electron quasiparticle density. The electronic excitations
at zero temperature necessarily carry the Fermi energy, and hence the 
wavefunction describing the excitations is of the Fermi wavelength. 
The electronic charge distrbution is the square modulus of the
wavefunction and hence takes on twice the periodicty or double the
wavelength $\lambda_\text{{Friedel}}=\lambda_\text{{Fermi}}/2=
{friedel_wavelength:.3}$. \\
''')

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
    plt.savefig(output+'.pdf', bbox_inches = "tight")

    with open(output+'.txt', 'w') as f:
        f.write(rf'''The local density of states of a spinless square lattice tight-binding
model at $\mu/t={mu:.2f}$ with ${nx}\times{ny}$
sites, a single orbital with an impurity at the centre with coupling
strength $V/t={V:.2f}$. The impurity gives rise to Fridel's eponymous
waves in the electron quasiparticle density. The electronic excitations
at zero temperature necessarily carry the Fermi energy, and hence the 
wavefunction describing the excitations is of the Fermi wavelength. 
The electronic charge distrbution is the square modulus of the
wavefunction and hence takes on twice the periodicty or double the
wavelength $\lambda_\text{{Friedel}}=\lambda_\text{{Fermi}}/2=
{friedel_wavelength:.3}$. \\
''')

def real_space(greens_function):
    FIGNAME='ldos'

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, axes=['resolved','integrated'],omega_min=-2,omega_max=2,vmin=0,vmax=80)

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH/2) 
    plt.tight_layout()
    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight")

    with open(output+'.txt', 'w') as f:
        f.write(rf'''The local density of states of a spinless square lattice tight-binding
model at $\mu/t={mu:.2f}$ with ${nx}\times{ny}$
sites, a single orbital with an impurity at the centre with coupling
strength $V/t={V:.2f}$. The impurity gives rise to Fridel's eponymous
waves in the electron quasiparticle density. The electronic excitations
at zero temperature necessarily carry the Fermi energy, and hence the 
wavefunction describing the excitations is of the Fermi wavelength. 
The electronic charge distrbution is the square modulus of the
wavefunction and hence takes on twice the periodicty or double the
wavelength $\lambda_\text{{Friedel}}=\lambda_\text{{Fermi}}/2=
{friedel_wavelength:.3}$. \\
''')

def k_space(greens_function):
    FIGNAME='k_space'

    fig, ax = plt.subplots()

    ax = greens_function.plot_spectrum(ax, axes=['integrated','resolved'],omega_min=-2,omega_max=2,vmin=0,vmax=80)

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH/2) 
    plt.tight_layout()
    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight")

    with open(output+'.txt', 'w') as f:
        f.write(rf'''The local density of states of a spinless square lattice tight-binding
model at $\mu/t={mu:.2f}$ with ${nx}\times{ny}$
sites, a single orbital with an impurity at the centre with coupling
strength $V/t={V:.2f}$. The impurity gives rise to Fridel's eponymous
waves in the electron quasiparticle density. The electronic excitations
at zero temperature necessarily carry the Fermi energy, and hence the 
wavefunction describing the excitations is of the Fermi wavelength. 
The electronic charge distrbution is the square modulus of the
wavefunction and hence takes on twice the periodicty or double the
wavelength $\lambda_\text{{Friedel}}=\lambda_\text{{Fermi}}/2=
{friedel_wavelength:.3}$. \\
''')

def majorana_fermi_arc(greens_function):
    FIGNAME='majorana_fermi_arc'

    fig, ax = plt.subplots()
    
    majorana_energy=0
    Bogoliubov_Fermi_arc_energy=0
    majorana_ky=2.12
    Bogoliubov_Fermi_arc_ky=0.5
    ax = greens_function.plot_spectrum(ax, energy=majorana_energy, axes=['resolved',majorana_ky], omega_min=0,omega_max='default',vmin=0,vmax=6,label='Majorana')
    ax = greens_function.plot_spectrum(ax, energy=Bogoliubov_Fermi_arc_energy, axes=['resolved',Bogoliubov_Fermi_arc_ky], omega_min=0,omega_max='default',vmin=0,vmax=6,label='Bogoloiubov-Fermi arc')
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
    plt.savefig(output+'.pdf', bbox_inches = "tight")

    # with open(output+'.txt', 'w') as f:
    #     f.write(rf'''The local density of states of a spinless square lattice tight-binding
# model at $\mu/t={mu:.2f}$ with ${nx}\times{ny}$
# sites, a single orbital with an impurity at the centre with coupling
# strength $V/t={V:.2f}$. The impurity gives rise to Fridel's eponymous
# waves in the electron quasiparticle density. The electronic excitations
# at zero temperature necessarily carry the Fermi energy, and hence the 
# wavefunction describing the excitations is of the Fermi wavelength. 
# The electronic charge distrbution is the square modulus of the
# wavefunction and hence takes on twice the periodicty or double the
# wavelength $\lambda_\text{{Friedel}}=\lambda_\text{{Fermi}}/2=
# {friedel_wavelength:.3}$. \\
# ''')
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
V=20
n_cells=41

# greens_function_xy, greens_function_xq, greens_function_kq, bdg = main()

[greens_function_xy, greens_function_xq, greens_function_kq, bdg] = np.load(DATA, allow_pickle=True)

# plot_itertions(bdg)
# unit_cell(bdg)
# ldos_each_atom(greens_function_xy)
# real_space(greens_function_xy)
# k_space(greens_function_kq)
majorana_fermi_arc(greens_function_xq)
