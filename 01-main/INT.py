from lib import *

filename=sys.argv[0].split('.')[0]
DATA_Uv=os.path.join(DATA,filename+'_Uv')
DATA_Delta=os.path.join(DATA,filename+'_Delta')
DATA=os.path.join(DATA,filename)
FIG=os.path.join(FIG,filename)
for directory in [DATA_Uv,FIG]:
    if not os.path.exists(directory):
        os.makedirs(directory)

def model(mu,Uv,s):
    """Creates a model for the phase diagram.
Uve two args: independent variables x and z"""
    A=Atom([0,0],'A')
    B=Atom([0.5,0],'B')
    A.add_orbital('s')
    B.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    bdg=BogoliubovdeGennes(lattice_vectors,'2D spinless p-wave theory')
    bdg.add_atom(A)
    bdg.add_atom(B)
    bdg.n_spins=1
    
    bdg.cut(n_cells, axes=0, glue_edgs=True)
    bdg.cut(n_cells, axes=1, glue_edgs=True)

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

    # bdg.add_impurities(V,[0,0])

    bdg.set_hartree(rho*1.1,atom='A')
    bdg.set_hartree(rho*0.9,atom='B')
    bdg.set_fock(phi,atom_i='A',atom_f='B')
    bdg.set_gorkov(chi,atom_i='A',atom_f='B')

    bdg.set_hubbard_u(-Uv,atom_i='A',atom_f='B',hop_vector=[0,0])

    bdg.record_hartree(location=[0,0],atom='A',_print=False)
    bdg.record_hartree(location=[0,0],atom='B',_print=False)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B', orbital_i='s', orbital_f='s', spin_i=0, spin_f=0,_print=False)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B', orbital_i='s', orbital_f='s', spin_i=0, spin_f=0,_print=False)

    return bdg

def model_Uv(Uv,mu):
    return model(mu,Uv,s)

def model_s(mu,s):
    return model(mu,Uv,s)

def dependent_variables(bdg):
    """The dependent variables to be extracted from the bdg model"""

    hartree=bdg.hartree()[0,0]
    fock=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])[0,0]
    gorkov=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])[0,0]

    return [hartree, fock, gorkov]

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
    ax1.plot(bdg._hartree_iterations[0],c='g',marker=markers[1],markersize=s,label=r'$\phi_{A}$', linestyle=linestyle)
    ax1.plot(bdg._hartree_iterations[1],c='g',marker=markers[1],markersize=s,label=r'$\phi_{B}$', linestyle=linestyle)
    ax1.plot(bdg._fock_iterations[0],c='k',marker=markers[2],markersize=s,label=r'$\Phi_{AB}$', linestyle=linestyle)
    ax1.plot(bdg._gorkov_iterations[0],c='cyan',marker=markers[3],markersize=s,label=r'$\Delta_{AB}$', linestyle=linestyle)
    ax1.legend()

    ax2.plot(bdg.free_energy,c='r',marker=markers[4],markersize=s, linestyle=linestyle)
    # plt.plot(np.real(bdg._gorkov_iterations[1]))
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    ax1.set_ylabel(r'Amplitude of fields')
    ax2.set_ylabel(r'Free energy')
    fig.suptitle(r'2D spinless p-wave theory')
    fig.supxlabel(r'Interaction $U$')
    
    plt.tight_layout()

    FIGNAME='2D spinless p-wave_renormalisation'

    output=os.path.join(FIG,FIGNAME)

    plt.savefig(output+'.pdf', bbox_inches = "tight")

    with open(output+'.txt', 'w') as f:
        f.write(rf'''Renormalisation of the Hartree, Fock and Gorkov fields on a 
        {n_cells}x{n_cells} Hubbard model with onsite repulsion $U={Uv}$.
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
    ax1.plot(bdg._fock_iterations[0],c='k',marker=markers[2],markersize=s,label=r'$\Phi_{AB}$')
    ax1.plot(bdg._gorkov_iterations[0],c='cyan',marker=markers[3],markersize=s,label=r'$\Delta_{AB}$')
    ax1.legend()

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('Free energy', color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.plot(bdg.free_energy,c='r',marker=markers[4],markersize=s,label=f'Free energy')

    ax3.plot(bdg._hartree_iterations[0],c='b',marker=markers[0],markersize=s,label=r'$\phi_{A}$')
    ax3.plot(bdg._hartree_iterations[1],c='g',marker=markers[1],markersize=s,label=r'$\phi_{B}$')
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

def plot_phase_diagram_Uv(phase_diagram_Uv):
    fig, [ax1,ax2] = plt.subplots(2,1,sharex=True)
    ax1 = phase_diagram_Uv.plot_phase_diagram(ax1,field_index=2)
    ax2 = phase_diagram_Uv.plot_phase_diagram(ax2,field_index=0)

    ax1.set_ylabel(r'$\Delta_{AB}$')
    ax2.set_ylabel('Free energy')
    ax2.set_xlabel(r'$U_v$')

    ax1.legend()
    ax1.get_legend().set_title(r'$\mu$')

    plt.subplots_adjust(hspace=.1)

    FIGNAME='Self-consistent_convergence_mu_Uv'

    fig.set_size_inches(w=LATEX_WIDTH, h=1.2*LATEX_WIDTH) 

    output=os.path.join(FIG,FIGNAME)
    plt.savefig(output+'.pdf', bbox_inches = "tight")

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
    fig, ax = plt.subplots(1,1)
    ax=phase_diagram_Delta.plot_phase_diagram(ax)
    ax.legend()

    FIGNAME='Self-consistent_convergence_nonzero_Delta'

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
mu=0.23
s=0
t=1
v=1.2
w=0.6
td=0.9
Uv=3.3
rho=-0.6
rho_shift=0.5
phi=0.6427
chi=1.828633
n_cells=11

# Phase diagram Uv, mu
rho_shift=3.443
rho=9.53378
phi=0.6427
chi=12.828633
init_friction=0.95
iter_friction=0.95
init_max_iterations=50
iter_max_iterations=30
absolute_convergence_factor=0.0001

_print=False
Uvv=np.arange(0.01,5,0.5)[::-1]
muu=np.arange(-5.,5,1.1)[::-1]
MUU=[3.5]

phase_diagram_Uv=PhaseDiagram(model_Uv)
phase_diagram_Uv.directory=DATA_Uv
phase_diagram_Uv.filename=filename
phase_diagram_Uv.initial_name='initial_convergence_Uv'
phase_diagram_Uv.initial_title=f'INT mean-field theory\n$\mu={muu[0]:.2f},\, U={Uvv[0]:.2f},\, \,$'
phase_diagram_Uv.phase_diagram(Uvv,dependent_variables,muu,init_friction=init_friction,iter_friction=iter_friction,init_max_iterations=init_max_iterations,iter_max_iterations=iter_max_iterations,absolute_convergence_factor=absolute_convergence_factor,initial_renormalisation_plot_function=plot_initial_renormalisation)
plot_phase_diagram_Uv(phase_diagram_Uv)

# phi=0.6427
# chi=1.828633
# phase_diagram_Delta=PhaseDiagram(model_Uv)
# phase_diagram_Delta.directory=DATA_Delta
# phase_diagram_Delta.filename=filename
# phase_diagram_Delta.initial_name='initial_convergence_nonzero_Delta'
# phase_diagram_Delta.initial_title=f'Self-consistent convergence\n$\mu={muu[0]:.2f},\, U={Uvv[0]:.2f},\, \,$'+rf'$\text{{friction}}={init_friction}$'
# # phase_diagram_Delta.phase_diagram(Uvv,dependent_variables,muu,init_friction=init_friction,iter_friction=iter_friction,init_max_iterations=init_max_iterations,iter_max_iterations=iter_max_iterations,absolute_convergence_factor=absolute_convergence_factor,initial_renormalisation_plot_function=plot_initial_renormalisation)
# plot_phase_diagram_Delta(phase_diagram_Delta)
