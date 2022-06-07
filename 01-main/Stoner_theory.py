from lib import *

filename=sys.argv[0].split('.')[0]
DATA_Us=os.path.join(DATA,filename+'_Us')
DATA_Delta=os.path.join(DATA,filename+'_Delta')
DATA=os.path.join(DATA,filename)
FIG=os.path.join(FIG,filename)
for directory in [DATA_Us,DATA_Delta,FIG]:
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

    bdg.cut(n_cells, axes=0, glue_edgs=True)
    bdg.cut(n_cells, axes=1, glue_edgs=True)
    
    bdg.set_onsite(-mu)

    bdg.set_hopping(-t,hop_vector=[1,0],label='$t$')
    bdg.set_hopping(-t,hop_vector=[0,1],label='$t$')

    bdg.set_hartree(rho+rho_shift,spin='up')
    bdg.set_hartree(rho-rho_shift,spin='dn')
    bdg.set_fock(phi,spin_i='up',spin_f='dn')
    bdg.set_gorkov(chi,spin_i='up',spin_f='dn')

    bdg.set_hubbard_u(Us,spin_i='up',spin_f='dn')

    bdg.record_hartree(location=[0,0], spin='up', _print=_print)
    bdg.record_hartree(location=[0,0], spin='dn', _print=_print)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='dn',_print=_print)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='dn',_print=_print)

    return bdg

def model_Us(Us,mu):
    return model(mu,Us,rho_shift)

def model_rho_shift(mu,rho_shift):
    return model(mu,Us,rho_shift)

def dependent_variables(bdg):
    """The dependent variables to be extracted from the bdg model"""

    hartree_up=bdg.hartree(spin='up')[0,0]
    hartree_dn=bdg.hartree(spin='dn')[0,0]
    fock=bdg.gorkov(spin_i='up', spin_f='dn', hop_vector=[0,0])[0,0]
    gorkov=bdg.gorkov(spin_i='up', spin_f='dn', hop_vector=[0,0])[0,0]

    return [hartree_up, hartree_dn, fock, gorkov]

def process(*args):
    """A function which processes the bdg model, returning greens functions"""

    bdg = model(*args)

    bdg.self_consistent_calculation(friction=friction, max_iterations=400, absolute_convergence_factor=absolute_convergence_factor)

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

    fig, [ax1, ax2] = plt.subplots(2,1,sharex='col')


    linestyle='dashed'

    # color = 'tab:black'
    ax1.tick_params(axis='y')
    magnetism=bdg._hartree_iterations[0]-bdg._hartree_iterations[1]
    ax1.plot(magnetism,c='g',marker=markers[1],markersize=s,label=r'$\hat{M}$', linestyle=linestyle)
    ax1.plot(bdg._fock_iterations[0],c='k',marker=markers[2],markersize=s,label=r'$\Phi_{\uparrow\downarrow}$', linestyle=linestyle)
    ax1.plot(bdg._gorkov_iterations[0],c='cyan',marker=markers[3],markersize=s,label=r'$\Delta_{\uparrow\downarrow}$', linestyle=linestyle)
    ax1.legend()

    ax2.plot(bdg.free_energy,c='r',marker=markers[4],markersize=s, linestyle=linestyle)
    # plt.plot(np.real(bdg._gorkov_iterations[1]))
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    ax1.set_ylabel(r'Amplitude of fields')
    ax2.set_ylabel(r'Free energy')
    fig.suptitle(r'Stoner theory')
    fig.supxlabel(r'Interaction $U$')
    
    plt.tight_layout()

    FIGNAME='Stoner_renormalisation'

    output=os.path.join(FIG,FIGNAME)

    plt.savefig(output+'.pdf', bbox_inches = "tight")

    with open(output+'.txt', 'w') as f:
        f.write(rf'''Renormalisation of the Hartree, Fock and Gorkov fields on a 
        {n_cells}x{n_cells} Hubbard model with onsite repulsion $U={Us}$.
attractive intracell Hubbard U, pairing fermions on sites A and B. 
The fields are calculated self-consistently, with absolute convergence factor 
${absolute_convergence_factor}$ and friction
${friction}$.
''')

def plot_initial_renormalisation(friction,max_iterations,absolute_convergence_factor,dataname,filename,title):

    markers=['o','+','^','x','.']
    s=3

    bdg = np.load(dataname+'.npz', allow_pickle=True)

    fig, [ax1, ax3] = plt.subplots(2,1,sharex='col')

    # color = 'tab:black'
    ax1.tick_params(axis='y')
    ax1.plot(bdg._fock_iterations[0],c='k',marker=markers[2],markersize=s,label=r'$\Phi_{\uparrow\downarrow}$')
    ax1.plot(bdg._gorkov_iterations[0],c='cyan',marker=markers[3],markersize=s,label=r'$\Delta_{\uparrow\downarrow}$')
    ax1.legend()

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('Free energy', color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.plot(bdg.free_energy,c='r',marker=markers[4],markersize=s,label=f'Free energy')

    ax3.plot(bdg._hartree_iterations[0],c='b',marker=markers[0],markersize=s,label=r'$\phi_{\uparrow}$')
    ax3.plot(bdg._hartree_iterations[1],c='g',marker=markers[1],markersize=s,label=r'$\phi_{\downarrow}$')
    ax3.legend()
    # plt.plot(np.real(bdg._gorkov_iterations[1]))
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    xlabel=r'Iterations'
    ylabel=r'Amplitude of fields'
    fig.suptitle(title)
    fig.supxlabel(xlabel)
    fig.supylabel(ylabel)
    plt.tight_layout()

    FIGNAME=filename

    output=os.path.join(FIG,FIGNAME)

    plt.savefig(output+'.pdf', bbox_inches = "tight")

    plt.close()

    with open(output+'.txt', 'w') as f:
        f.write(rf'''Renormalisation of the Hartree, Fock and Gorkov fields on a 
        2-dimensional spinless Weyl-SSH lattice with 
attractive intracell Hubbard U, pairing fermions on sites A and B. 
The fields are calculated self-consistently, with absolute convergence factor 
${absolute_convergence_factor}$
and an initial friction ${init_friction}$ and 
friction ${iter_friction}$ over subsequent iterations.
The lattice is ${n_cells}x{n_cells}$ cells squared.
''')

def plot_phase_diagram_Us(phase_diagram_Us):
    fig, axs = plt.subplots(5,1,sharex=True)
    axs[0] = phase_diagram_Us.plot_phase_diagram(axs[0],field_index=1,minus_field_index=2)
    axs[1] = phase_diagram_Us.plot_phase_diagram(axs[1],field_index=1,plus_field_index=2)
    axs[2] = phase_diagram_Us.plot_phase_diagram(axs[2],field_index=3)
    axs[3] = phase_diagram_Us.plot_phase_diagram(axs[3],field_index=4)
    axs[4] = phase_diagram_Us.plot_phase_diagram(axs[4],field_index=0)

    axs[0].set_ylabel(r'$\langle\hat{M}\rangle$')
    axs[1].set_ylabel(r'$\sum_\sigma\phi_\sigma$')
    axs[2].set_ylabel(r'$\Phi_{\uparrow\downarrow}$')
    axs[3].set_ylabel(r'$\Delta_{\uparrow\downarrow}$')
    axs[4].set_ylabel('Free energy')
    axs[4].set_xlabel(r'$U_s$')

    axs[0].legend(loc='upper left', bbox_to_anchor=(1.05,1))
    axs[0].get_legend().set_title(r'$\mu$')

    # plt.subplots_adjust(hspace=.1)

    FIGNAME='Self-consistent_convergence_mu_Us'

    fig.set_size_inches(w=LATEX_WIDTH, h=1.2*LATEX_WIDTH) 

    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight")
    plt.close()

    with open(output+'.txt', 'w') as f:
        f.write(rf'''Self-consistent Hartree, Fock and Gorkov fields on an 
${n_cells}x{n_cells}$ lattice as a function of Coulomb repulsion U for 
various chemical potentials $\mu$.
The absolute convergence factor is
${absolute_convergence_factor}$
and an initial friction ${init_friction}$ and 
friction ${iter_friction}$ over subsequent iterations.
''')

def plot_phase_diagram_Delta(phase_diagram_Delta):
    fig, axs = plt.subplots(5,1,sharex=True)
    axs[0] = phase_diagram_Delta.plot_phase_diagram(axs[0],field_index=1,minus_field_index=2)
    axs[1] = phase_diagram_Delta.plot_phase_diagram(axs[1],field_index=1,plus_field_index=2)
    axs[2] = phase_diagram_Delta.plot_phase_diagram(axs[2],field_index=3)
    axs[3] = phase_diagram_Delta.plot_phase_diagram(axs[3],field_index=4)
    axs[-1] = phase_diagram_Delta.plot_phase_diagram(axs[-1],field_index=0)

    axs[0].set_ylabel(r'$\langle\hat{M}\rangle$')
    axs[1].set_ylabel(r'$\sum_\sigma\phi_\sigma$')
    axs[2].set_ylabel(r'$\Phi_{\uparrow\downarrow}$')
    axs[3].set_ylabel(r'$\Delta_{\uparrow\downarrow}$')
    axs[-1].set_ylabel('Free energy')
    axs[-1].set_xlabel(r'$U_s$')

    axs[0].legend(loc='upper left', bbox_to_anchor=(1.05,1))
    axs[0].get_legend().set_title(r'$\mu$')

    # plt.subplots_adjust(hspace=1)

    FIGNAME='Self-consistent_convergence_Delta'

    fig.set_size_inches(w=LATEX_WIDTH, h=1.2*LATEX_WIDTH) 

    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight")
    plt.close()

    with open(output+'.txt', 'w') as f:
        f.write(rf'''Self-consistent Hartree, Fock and Gorkov fields on an 
${n_cells}x{n_cells}$ lattice as a function of Coulomb repulsion U for 
various chemical potentials $\mu$.
The absolute convergence factor is
${absolute_convergence_factor}$
and an initial friction ${init_friction}$ and 
friction ${iter_friction}$ over subsequent iterations.
''')

def plot_phase_diagram(phase_diagram_Us,phase_diagram_Delta):
    fig, axs = plt.subplots(5,2,sharex=True,sharey='row')
    axs[0,0] = phase_diagram_Us.plot_phase_diagram(axs[0,0],field_index=1,minus_field_index=2)
    axs[1,0] = phase_diagram_Us.plot_phase_diagram(axs[1,0],field_index=1,plus_field_index=2)
    axs[2,0] = phase_diagram_Us.plot_phase_diagram(axs[2,0],field_index=3)
    axs[3,0] = phase_diagram_Us.plot_phase_diagram(axs[3,0],field_index=4)
    axs[4,0] = phase_diagram_Us.plot_phase_diagram(axs[4,0],field_index=0)

    axs[0,0].set_ylabel(r'$\langle\hat{M}\rangle$')
    axs[1,0].set_ylabel(r'$\sum_\sigma\phi_\sigma$')
    axs[2,0].set_ylabel(r'$\Phi_{\uparrow\downarrow}$')
    axs[3,0].set_ylabel(r'$\Delta_{\uparrow\downarrow}$')
    axs[4,0].set_ylabel('Free energy')
    # fig.suptitle(greens_function.title)
    fig.supxlabel(r'$U_s$')

    axs[0,1] = phase_diagram_Delta.plot_phase_diagram(axs[0,1],field_index=1,minus_field_index=2)
    axs[1,1] = phase_diagram_Delta.plot_phase_diagram(axs[1,1],field_index=1,plus_field_index=2)
    axs[2,1] = phase_diagram_Delta.plot_phase_diagram(axs[2,1],field_index=2)
    axs[3,1] = phase_diagram_Delta.plot_phase_diagram(axs[3,1],field_index=3)
    axs[-1,1] = phase_diagram_Delta.plot_phase_diagram(axs[-1,1],field_index=0)

    axs[0,1].legend(loc='upper left', bbox_to_anchor=(1.05,1))
    axs[0,1].get_legend().set_title(r'$\mu$')

    axs[0,0].set_title('Initial magnetism')
    axs[0,1].set_title('Initial BCS pairing')

    # plt.subplots_adjust(hspace=.1)

    FIGNAME='Self-consistent_convergence'

    fig.set_size_inches(w=LATEX_WIDTH, h=1.2*LATEX_WIDTH) 

    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight")
    plt.close()

    with open(output+'.txt', 'w') as f:
        f.write(rf'''Self-consistent Hartree, Fock and Gorkov fields on an 
${n_cells}x{n_cells}$ lattice as a function of Coulomb repulsion U for 
various chemical potentials $\mu$.
The absolute convergence factor is
${absolute_convergence_factor}$
and an initial friction ${init_friction}$ and 
friction ${iter_friction}$ over subsequent iterations.
''')

#############################################################################
################################# Main ######################################
#############################################################################
mu=0.23
t=1
Us=3.3
rho=-0.6
rho_shift=0.5
phi=0.6427
chi=1.828633
n_cells=43

friction=0.95
absolute_convergence_factor=0.00001

bulk_calculation=True

_print=False
bdg = model(mu,Us,rho_shift)

# bdg.self_consistent_calculation(friction=friction, max_iterations=100, absolute_convergence_factor=absolute_convergence_factor)

# hartree_up=bdg.hartree(spin='up')
# hartree_dn=bdg.hartree(spin='dn')
# fock=bdg.fock(spin_i='up',spin_f='dn')
# gorkov=bdg.gorkov(spin_i='up',spin_f='dn')
# M=(hartree_up-hartree_dn)[0,0]
# plot_iterations(bdg)
# plt.show()
# exit()

# Phase diagram Us, mu
rho_shift=3.443
rho=-2.53378
phi=0 #0.6427
chi=0 #1.828633
init_friction=0
iter_friction=0
init_max_iterations=50
iter_max_iterations=20
absolute_convergence_factor=0.0001

plot_initial_renormalisation=None

_print=False
Uss=np.arange(-10,10,0.51)[::-1]
muu=np.arange(-5,5.1,1.0)[::-1]
# muu=np.arange(-5,-3.1,1.0)[::-1]

phase_diagram_Us=PhaseDiagram(model_Us)
phase_diagram_Us.directory=DATA_Us
phase_diagram_Us.filename=filename
phase_diagram_Us.initial_name='initial_convergence_Us'
phase_diagram_Us.initial_title=f'Stoner mean-field theory\n$\mu={muu[0]:.2f},\, U={Uss[0]:.2f},\, \,$'
phase_diagram_Us.phase_diagram(Uss,dependent_variables,muu,init_friction=init_friction,iter_friction=iter_friction,init_max_iterations=init_max_iterations,iter_max_iterations=iter_max_iterations,absolute_convergence_factor=absolute_convergence_factor,initial_renormalisation_plot_function=plot_initial_renormalisation,plot_phase_diagram_function=plot_phase_diagram_Us)
# plot_phase_diagram_Us(phase_diagram_Us)

iter_max_iterations=50
rho=-6.53378
rho_shift=0#.443
phi=0
chi=5.828633
phase_diagram_Delta=PhaseDiagram(model_Us)
phase_diagram_Delta.directory=DATA_Delta
phase_diagram_Delta.filename=filename
phase_diagram_Delta.initial_name='initial_convergence_Delta'
phase_diagram_Delta.initial_title=f'Stoner theory with nonzero Gorkov pairing\n$\mu={muu[0]:.2f},\, U={Uss[0]:.2f},\, \,$'
phase_diagram_Delta.phase_diagram(Uss,dependent_variables,muu,init_friction=init_friction,iter_friction=iter_friction,init_max_iterations=init_max_iterations,iter_max_iterations=iter_max_iterations,absolute_convergence_factor=absolute_convergence_factor,initial_renormalisation_plot_function=plot_initial_renormalisation,plot_phase_diagram_function=plot_phase_diagram_Delta)
# plot_phase_diagram_Delta(phase_diagram_Delta)

plot_phase_diagram(phase_diagram_Us,phase_diagram_Delta)
