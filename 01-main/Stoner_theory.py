from lib import *

filename=sys.argv[0].split('.')[0]
DATA_Uv=os.path.join(DATA,filename+'_Uv')
DATA_s=os.path.join(DATA,filename+'_s')
DATA=os.path.join(DATA,filename)
FIG=os.path.join(FIG,filename)
for directory in [DATA_Uv,DATA_s,FIG]:
    if not os.path.exists(directory):
        os.makedirs(directory)

def model(mu,Us,rho_shift):
    """Creates a model for the phase diagram.
Use two args: independent variables x and z"""
    A=Atom([0,0],'A')
    A.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    bdg=BogoliubovdeGennes(lattice_vectors,'2D-Weyl-SSH')
    bdg.add_atom(A)
    bdg.n_spins=2

    # bdg.set_kpts([n_cells,n_cells])

    bdg.cut(n_cells, axes=0, glue_edgs=True)
    bdg.cut(n_cells, axes=1, glue_edgs=True)
    bdg.set_onsite(-mu)

    bdg.set_hopping(-t,hop_vector=[1,0],label='$t$')
    bdg.set_hopping(-t,hop_vector=[0,1],label='$t$')

    bdg.set_hartree(rho+rho_shift,spin='up')
    bdg.set_hartree(rho-rho_shift,spin='dn')
    bdg.set_fock(phi,spin_i='up',spin_f='dn')
    bdg.set_gorkov(chi,spin_i='up',spin_f='dn')

    bdg.set_hubbard_u(-Us,spin_i='up',spin_f='dn')
    
    _print=True
    bdg.record_hartree(location=[0,0], spin='up', _print=_print)
    bdg.record_hartree(location=[0,0], spin='dn', _print=_print)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='dn',_print=_print)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='dn',_print=_print)

    return bdg

def model_Us(mu,Us):
    return model(mu,Us,rho_shift)

def model_rho_shift(mu,rho_shift):
    return model(mu,Us,rho_shift)

def dependent_variables(bdg):
    """The dependent variables to be extracted from the bdg model"""

    hartree_up=bdg.hartree(spin='up')[0,0]
    hartree_dn=bdg.hartree(spin='dn')[0,0]
    fock=bdg.gorkov(spin_i='up', spin_f='dn', hop_vector=[0,0])[0,0]
    gorkov=bdg.gorkov(spin_i='up', spin_f='dn', hop_vector=[0,0])[0,0]

    return [hartree_up, hartree_dn, fock_v, gorkov_v]

def process(*args):
    """A function which processes the bdg model, returning greens functions"""

    bdg = model(*args)

    bdg.self_consistent_calculation(friction=0.9, max_iterations=400, absolute_convergence_factor=0.00001)

    energy_interval=np.linspace(-4,4,601)
    resolution=0.05

    # greens_function_xy=GreensFunction(bdg,energy_interval,resolution, k_axes=None)

    # greens_function_xq=GreensFunction(bdg,energy_interval,resolution, k_axes=[1])

    greens_function_kq=GreensFunction(bdg,energy_interval,resolution, k_axes=[0,1])

    del bdg.eigenvectors
    del bdg.eigenvalues

    with open(DATA+'.npz', 'wb') as f:
        cPickle.dump([greens_function_xy, greens_function_xq, greens_function_kq, bdg], f)
    return greens_function_kq, bdg
    return greens_function_xy, greens_function_xq, greens_function_kq, bdg

# Plotting:

def plot_iterations(bdg):

    markers=['o','+','^','x','.']
    s=3

    # hartree_A=bdg.hartree(atom='A')[0,0]
    # hartree_B=bdg.hartree(atom='B')[0,0]
    # fock_v=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])[0,0]
    # gorkov_v=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])[0,0]
    # free_energy=bdg.free_energy
    
    # exit()
    # gorkov_w=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[-1,0])

    fig, [ax1, ax3] = plt.subplots(2,1,sharex='col')

    # color = 'tab:black'
    ax1.tick_params(axis='y')
    ax1.plot(bdg._fock_iterations[0],c='k',marker=markers[2],markersize=s,label=f'$\Phi$')
    ax1.plot(bdg._gorkov_iterations[0],c='cyan',marker=markers[3],markersize=s,label=f'$\Delta$')
    ax1.legend()

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('Free energy', color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.plot(bdg.free_energy,c='r',marker=markers[4],markersize=s,label=f'Free energy')
    # ax2.plot(bdg.Eg,c='b',marker=markers[4],markersize=s,label=f'Eg')
    # ax2.plot(bdg.V,c='g',marker=markers[4],markersize=s,label=f'V')
    # ax2.plot(bdg.V_mf,c='c',marker=markers[4],markersize=s,label=r'V_{mf}')

    ax3.plot(bdg._hartree_iterations[0],c='b',marker=markers[0],markersize=s,label=r'$\phi_\uparrow$')
    ax3.plot(bdg._hartree_iterations[1],c='g',marker=markers[1],markersize=s,label=r'$\phi_\downarrow$')
    ax3.legend()
    # plt.plot(np.real(bdg._gorkov_iterations[1]))
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    xlabel=r'Iterations'
    ylabel=r'Amplitude of fields'
    title=r'Stoner theory'
    fig.suptitle(title)
    fig.supxlabel(xlabel)
    fig.supylabel(ylabel)
    plt.tight_layout()
    
    # print('Hartree A:')
    # print(bdg._hartree_iterations[0,-1])
    # print('Hartree B:')
    # print(bdg._hartree_iterations[1,-1])
    # print('Fock v:')
    # print(bdg._fock_iterations[0,-1])
    # print('Gorkov v:')
    # print(bdg._gorkov_iterations[0,-1])

    # FIGNAME=filename
    plt.show()

    # output=os.path.join(FIG,FIGNAME)

    # plt.savefig(output+'.pdf', bbox_inches = "tight")

    # with open(output+'.txt', 'w') as f:
    #     f.write(rf'''Renormalisation of the Hartree, Fock and Gorkov fields on a 
    #     2-dimensional spinless Weyl-SSH lattice with 
# attractive intracell Hubbard U, pairing fermions on sites A and B. 
# The fields are calculated self-consistently, with absolute convergence factor 
# ${absolute_convergence_factor}$
# and an initial friction ${init_friction}$ and 
# friction ${iter_friction}$ over subsequent iterations.
# The lattice is ${n_cells}x{n_cells}$ cells squared.
# ''')

def plot_phase_diagram_Uv(phase_diagram_Uv):
    fig, ax = plt.subplots(1,1)
    ax=phase_diagram_Uv.plot_phase_diagram(ax)
    ax.legend()

    FIGNAME='Self-consistent_convergence_mu_Uv'

    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight")

    with open(output+'.txt', 'w') as f:
        f.write(rf'''Self-consistent Hartree, Fock and Gorkov fields on an 
${n_cells}x{n_cells}$ 
spinless Weyl-SSH lattice as a function of chemical potential 
$\mu$ and multiorbital Hubbard attraction $U_v$.
The absolute convergence factor is
${absolute_convergence_factor}$
and an initial friction ${init_friction}$ and 
friction ${iter_friction}$ over subsequent iterations.
''')

def plot_phase_diagram_s(phase_diagram_s):
    fig, ax = plt.subplots(1,1)
    ax=phase_diagram_s.plot_phase_diagram(ax)
    ax.legend()

    FIGNAME='Self-consistent_convergence_mu_s'

    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight")

    with open(output+'.txt', 'w') as f:
        f.write(rf'''Self-consistent Hartree, Fock and Gorkov fields on an ${n_cells}x{n_cells}$ spinless Weyl-SSH lattice as a function of chemical potential $\mu={mu}$ and multiorbital rigid intracell shift $s$, with $U_v={Uv}$.
The absolute convergence factor is
${absolute_convergence_factor}$
and an initial friction ${init_friction}$ and 
friction ${iter_friction}$ over subsequent iterations.
''')

#############################################################################
################################# Main ######################################
#############################################################################
mu=4.8
t=1
Us=-5.8
rho=-4.15
rho_shift=1.15
phi=0.
chi=0#3.4
n_cells=11


bdg = model_Us(mu,Us)

bdg.self_consistent_calculation(friction=0.8, max_iterations=400, absolute_convergence_factor=0.0001)
# exit()

hartree_up=bdg.hartree(spin='up')[0,0]
hartree_dn=bdg.hartree(spin='dn')[0,0]
M=hartree_up-hartree_dn
print(M)
# exit()
# Convergence plot
plot_iterations(bdg)
exit()

# Phase diagram mu,Uv/s
rho_shift=0.4
rho=-2.13378
phi_v=0.06427
chi_v=0.828633
init_friction=0.6
iter_friction=0.7
init_max_iterations=400
iter_max_iterations=200
absolute_convergence_factor=0.0001

muu=np.arange(-4,4.1,0.1)[::-1]
Uvv=np.arange(1.1,6.5,1.1)[::-1]
phase_diagram_Uv=PhaseDiagram(model_Uv)
phase_diagram_Uv.directory=DATA_Uv
phase_diagram_Uv.filename=filename
phase_diagram_Uv.initial_name='initial_convergence_Uv'
phase_diagram_Uv.initial_title=f'Self-consistent convergence\n$U_v={Uvv[0]:.2f},\, \mu={muu[0]:.2f},\, s={s},\,$'+rf'$\text{{friction}}={init_friction}$'
phase_diagram_Uv.phase_diagram(muu,dependent_variables,Uvv,init_friction=init_friction,iter_friction=iter_friction,init_max_iterations=init_max_iterations,iter_max_iterations=iter_max_iterations,absolute_convergence_factor=absolute_convergence_factor,initial_renormalisation_plot_function=plot_initial_renormalisation)
plot_phase_diagram_Uv(phase_diagram_Uv)


Uv=3.6
ss=np.arange(0,2,0.44)[::-1]
phase_diagram_s=PhaseDiagram(model_s)
phase_diagram_s.directory=DATA_s
phase_diagram_s.initial_name='initial_convergence_s'
phase_diagram_s.filename=filename
phase_diagram_s.initial_title=f'Self-consistent convergence\n$U_v={Uv},\, \mu={muu[0]:.2f},\, s={ss[0]:.2f},\,$'+rf'$\text{{friction}}={init_friction}$'
phase_diagram_s.phase_diagram(muu,dependent_variables,ss,init_friction=init_friction,iter_friction=iter_friction,init_max_iterations=init_max_iterations,iter_max_iterations=iter_max_iterations,absolute_convergence_factor=absolute_convergence_factor,initial_renormalisation_plot_function=plot_initial_renormalisation)
plot_phase_diagram_s(phase_diagram_s)
