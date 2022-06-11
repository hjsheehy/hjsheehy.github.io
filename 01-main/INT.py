from lib import *

filename=sys.argv[0].split('.')[0]
DATA=os.path.join(DATA,filename)
FIG=os.path.join(FIG,filename)
for directory in [DATA,FIG]:
    if not os.path.exists(directory):
        os.makedirs(directory)

def model(mu,Us,Uv,s):
    """Creates a model for the phase diagram.
Uve two args: independent variables x and z"""
    A=Atom([0,0],'A')
    B=Atom([0.5,0],'B')
    A.add_orbital('s')
    B.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    bdg=BogoliubovdeGennes(lattice_vectors,'Non-unitary triplet theory')
    bdg.add_atom(A)
    bdg.add_atom(B)
    bdg.n_spins=2
    
    bdg.cut(n_cells, axes=0, glue_edgs=True)
    bdg.cut(n_cells, axes=1, glue_edgs=True)

    bdg.set_onsite(-mu+s,atom='A')
    bdg.set_onsite(-mu-s,atom='B')

    # bdg.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B',label='$v$')
    # bdg.set_hopping(-w,hop_vector=[-1,0],atom_i='A',atom_f='B',label='$w$')
    # bdg.set_hopping(-td,hop_vector=[0,1],atom_i='B',atom_f='A',label='$t_d$')
    # bdg.set_hopping(-td,hop_vector=[0,-1],atom_i='B',atom_f='A',label='$t_d$')
    # bdg.set_hopping(-td,hop_vector=[1,1],atom_i='B',atom_f='A',label='$t_d$')
    # bdg.set_hopping(-td,hop_vector=[1,-1],atom_i='B',atom_f='A',label='$t_d$')
    
    # impurity_wall = [[0,i] for i in range(n_cells)]
    # bdg.add_impurities(V,impurity_wall)

    # bdg.add_impurities(V,[0,0])

    bdg.set_hartree(rho+rho_shift,atom='A',spin='up')
    bdg.set_hartree(rho+rho_shift,atom='A',spin='dn')
    bdg.set_hartree(rho-rho_shift,atom='B',spin='up')
    bdg.set_hartree(rho-rho_shift,atom='B',spin='dn')
    bdg.set_fock(phi,atom_i='A',atom_f='B',spin_i='up',spin_f='up')
    bdg.set_fock(phi,atom_i='A',atom_f='B',spin_i='dn',spin_f='dn')
    bdg.set_fock(phi,atom_i='A',atom_f='A',spin_i='up',spin_f='dn')
    bdg.set_fock(phi,atom_i='B',atom_f='B',spin_i='up',spin_f='dn')
    bdg.set_gorkov(chi+chi_shift,atom_i='A',atom_f='B',spin_i='up',spin_f='up')
    bdg.set_gorkov(chi-chi_shift,atom_i='A',atom_f='B',spin_i='dn',spin_f='dn')
    bdg.set_gorkov(chi,atom_i='A',atom_f='A',spin_i='up',spin_f='dn')
    bdg.set_gorkov(chi,atom_i='B',atom_f='B',spin_i='up',spin_f='dn')
    
    bdg.set_hubbard_u(-Uv,atom_i='A',atom_f='B',spin_i='up',spin_f='up',hop_vector=[0,0])
    bdg.set_hubbard_u(-Uv,atom_i='A',atom_f='B',spin_i='dn',spin_f='dn',hop_vector=[0,0])
    bdg.set_hubbard_u(Us,atom_i='A',atom_f='A',spin_i='up',spin_f='dn',hop_vector=[0,0])
    bdg.set_hubbard_u(Us,atom_i='B',atom_f='B',spin_i='up',spin_f='dn',hop_vector=[0,0])
    
    _print=False
    bdg.record_hartree(location=[0,0],atom='A',spin='up',_print=_print)
    bdg.record_hartree(location=[0,0],atom='A',spin='dn',_print=_print)
    bdg.record_hartree(location=[0,0],atom='B',spin='up',_print=_print)
    bdg.record_hartree(location=[0,0],atom='B',spin='dn',_print=_print)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], atom_i='A',atom_f='B',spin_i='up',spin_f='up')
    bdg.record_fock(location_i=[0,0], location_f=[0,0], atom_i='A',atom_f='B',spin_i='dn',spin_f='dn')
    bdg.record_fock(location_i=[0,0], location_f=[0,0], atom_i='A',atom_f='A',spin_i='up',spin_f='dn')
    bdg.record_fock(location_i=[0,0], location_f=[0,0], atom_i='B',atom_f='B',spin_i='up',spin_f='dn')
    # _print=True
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B', spin_i=0, spin_f=0,_print=_print)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B', spin_i=1, spin_f=1,_print=_print)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='A', spin_i=0, spin_f=1,_print=_print)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='B', atom_f='B', spin_i=0, spin_f=1,_print=_print)

    return bdg

def model_Uv(mu,Uv):
    return model(mu,Us,Uv,s)

def model_s(mu,s):
    return model(mu,Us,Uv,s)

def dependent_variables(bdg):
    """The dependent variables to be extracted from the bdg model"""

    hartree_A_up=bdg.hartree(atom='A',spin='up')[0,0]
    hartree_A_dn=bdg.hartree(atom='A',spin='dn')[0,0]
    hartree_B_up=bdg.hartree(atom='B',spin='up')[0,0]
    hartree_B_dn=bdg.hartree(atom='B',spin='dn')[0,0]
    fock_A_B_up_up=bdg.fock(atom_i='A', atom_f='B', spin_i='up',spin_f='up', hop_vector=[0,0])[0,0]
    fock_A_B_dn_dn=bdg.fock(atom_i='A', atom_f='B', spin_i='dn',spin_f='dn', hop_vector=[0,0])[0,0]
    fock_A_A_up_dn=bdg.fock(atom_i='A', atom_f='A', spin_i='up',spin_f='dn', hop_vector=[0,0])[0,0]
    fock_B_B_up_dn=bdg.fock(atom_i='B', atom_f='B', spin_i='up',spin_f='dn', hop_vector=[0,0])[0,0]
    gorkov_A_B_up_up=bdg.gorkov(atom_i='A', atom_f='B', spin_i='up',spin_f='up', hop_vector=[0,0])[0,0]
    gorkov_A_B_dn_dn=bdg.gorkov(atom_i='A', atom_f='B', spin_i='dn',spin_f='dn', hop_vector=[0,0])[0,0]
    gorkov_A_A_up_dn=bdg.gorkov(atom_i='A', atom_f='A', spin_i='up',spin_f='dn', hop_vector=[0,0])[0,0]
    gorkov_B_B_up_dn=bdg.gorkov(atom_i='B', atom_f='B', spin_i='up',spin_f='dn', hop_vector=[0,0])[0,0]
    return [hartree_A_up, hartree_A_dn, hartree_B_up, hartree_B_dn, 
            fock_A_B_up_up, fock_A_B_dn_dn, fock_A_A_up_dn, fock_B_B_up_dn, 
            gorkov_A_B_up_up, gorkov_A_B_dn_dn, gorkov_A_A_up_dn, gorkov_B_B_up_dn]

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

    hartree_A_up=bdg.hartree(atom='A',spin='up')[0,0]
    hartree_A_dn=bdg.hartree(atom='A',spin='dn')[0,0]
    hartree_B_up=bdg.hartree(atom='B',spin='up')[0,0]
    hartree_B_dn=bdg.hartree(atom='B',spin='dn')[0,0]
    fock_A_B_up_up=bdg.fock(atom_i='A', atom_f='B', spin_i='up',spin_f='up', hop_vector=[0,0])[0,0]
    fock_A_B_dn_dn=bdg.fock(atom_i='A', atom_f='B', spin_i='dn',spin_f='dn', hop_vector=[0,0])[0,0]
    fock_A_A_up_dn=bdg.fock(atom_i='A', atom_f='A', spin_i='up',spin_f='dn', hop_vector=[0,0])[0,0]
    fock_B_B_up_dn=bdg.fock(atom_i='B', atom_f='B', spin_i='up',spin_f='dn', hop_vector=[0,0])[0,0]
    gorkov_A_B_up_up=bdg.gorkov(atom_i='A', atom_f='B', spin_i='up',spin_f='up', hop_vector=[0,0])[0,0]
    gorkov_A_B_dn_dn=bdg.gorkov(atom_i='A', atom_f='B', spin_i='dn',spin_f='dn', hop_vector=[0,0])[0,0]
    gorkov_A_A_up_dn=bdg.gorkov(atom_i='A', atom_f='A', spin_i='up',spin_f='dn', hop_vector=[0,0])[0,0]
    gorkov_B_B_up_dn=bdg.gorkov(atom_i='B', atom_f='B', spin_i='up',spin_f='dn', hop_vector=[0,0])[0,0]
    free_energy=bdg.free_energy

    fig, [ax1, ax2, ax3] = plt.subplots(3,1,sharex='col')

    linestyle=['solid','dashed','dotted','dashdot']

    # color = 'tab:black'
    ax1.tick_params(axis='y')
    magnetism_A=bdg._hartree_iterations[0]-bdg._hartree_iterations[1]
    magnetism_B=bdg._hartree_iterations[2]-bdg._hartree_iterations[3]
    ax1.plot(magnetism_A,c='g',marker=markers[0],markersize=s,label=r'$\langle\hat{M}_A\rangle$', linestyle=linestyle[0])
    ax1.plot(magnetism_B,c='b',marker=markers[1],markersize=s,label=r'$\langle\hat{M}_B\rangle$', linestyle=linestyle[1])
    ax2.plot(bdg._fock_iterations[0],c='r',marker=markers[0],markersize=s,label=r'$\Phi^{AB}_{\uparrow\uparrow}$', linestyle=linestyle[0])
    ax2.plot(bdg._fock_iterations[1],c='g',marker=markers[1],markersize=s,label=r'$\Phi^{AB}_{\downarrow\downarrow}$', linestyle=linestyle[1])
    ax2.plot(bdg._fock_iterations[2],c='b',marker=markers[2],markersize=s,label=r'$\Phi^{AA}_{\uparrow\downarrow}$', linestyle=linestyle[2])
    ax2.plot(bdg._fock_iterations[3],c='k',marker=markers[3],markersize=s,label=r'$\Phi^{BB}_{\uparrow\downarrow}$', linestyle=linestyle[3])

    ax3.plot(bdg._gorkov_iterations[0],c='r',marker=markers[0],markersize=s,label=r'$\Delta^{AB}_{\uparrow\uparrow}$', linestyle=linestyle[0])
    ax3.plot(bdg._gorkov_iterations[1],c='g',marker=markers[1],markersize=s,label=r'$\Delta^{AB}_{\downarrow\downarrow}$', linestyle=linestyle[1])
    ax3.plot(bdg._gorkov_iterations[2],c='b',marker=markers[2],markersize=s,label=r'$\Delta^{AA}_{\uparrow\downarrow}$', linestyle=linestyle[2])
    ax3.plot(bdg._gorkov_iterations[3],c='k',marker=markers[3],markersize=s,label=r'$\Delta^{BB}_{\uparrow\downarrow}$', linestyle=linestyle[3])


    ax4 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax4.set_ylabel('Free energy', color=color)  # we already handled the x-label with ax1
    ax4.tick_params(axis='y', labelcolor=color)
    ax4.plot(bdg.free_energy,c='r',marker=markers[4],markersize=s,label=f'Free energy')

    ax1.legend()
    ax2.legend()
    ax3.legend()

    # plt.plot(np.real(bdg._gorkov_iterations[1]))
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    title=r'''2D Multiorbital non-unitary spin-triplet theory
'''+rf'$\mu={mu}, v={v}, w={w}, t_d={td}$'+r'''
'''+rf'$U_v={Uv}, U_s={Us}$'

    fig.suptitle(title)
    fig.supxlabel(r'Iterations')
    fig.supylabel(r'Amplitude of fields')
    
    plt.tight_layout()
    
    FIGNAME='INT_renormalisaiton'
    OUTPUT=os.path.join(FIG,FIGNAME)

    plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight")

    with open(OUTPUT+'.txt', 'w') as f:
        f.write(rf'''Renormalisation of the Hartree, Fock and Gorkov fields on a 
        {n_cells}x{n_cells} Hubbard model with onsite repulsion $U={Uv}$.
attractive intracell Hubbard U, pairing fermions on sites A and B. 
The fields are calculated self-consistently, with absolute convergence factor 
${absolute_convergence_factor}$ and friction
${friction}$.
''')

def plot_iterations_uni_vs_non(bdg_uni,bdg_non):

    markers=['o','+','^','x','.']
    s=3

    hartree_A_up=bdg_uni.hartree(atom='A',spin='up')[0,0]
    hartree_A_dn=bdg_uni.hartree(atom='A',spin='dn')[0,0]
    hartree_B_up=bdg_uni.hartree(atom='B',spin='up')[0,0]
    hartree_B_dn=bdg_uni.hartree(atom='B',spin='dn')[0,0]
    fock_A_B_up_up=bdg_uni.fock(atom_i='A', atom_f='B', spin_i='up',spin_f='up', hop_vector=[0,0])[0,0]
    fock_A_B_dn_dn=bdg_uni.fock(atom_i='A', atom_f='B', spin_i='dn',spin_f='dn', hop_vector=[0,0])[0,0]
    fock_A_A_up_dn=bdg_uni.fock(atom_i='A', atom_f='A', spin_i='up',spin_f='dn', hop_vector=[0,0])[0,0]
    fock_B_B_up_dn=bdg_uni.fock(atom_i='B', atom_f='B', spin_i='up',spin_f='dn', hop_vector=[0,0])[0,0]
    gorkov_A_B_up_up=bdg_uni.gorkov(atom_i='A', atom_f='B', spin_i='up',spin_f='up', hop_vector=[0,0])[0,0]
    gorkov_A_B_dn_dn=bdg_uni.gorkov(atom_i='A', atom_f='B', spin_i='dn',spin_f='dn', hop_vector=[0,0])[0,0]
    gorkov_A_A_up_dn=bdg_uni.gorkov(atom_i='A', atom_f='A', spin_i='up',spin_f='dn', hop_vector=[0,0])[0,0]
    gorkov_B_B_up_dn=bdg_uni.gorkov(atom_i='B', atom_f='B', spin_i='up',spin_f='dn', hop_vector=[0,0])[0,0]
    free_energy=bdg_uni.free_energy

    fig, axs = plt.subplots(4,2,sharex='col',sharey='row')

    linestyle=['solid','dashed','dotted','dashdot']

    # color = 'tab:black'
    axs[0,0].tick_params(axis='y')
    magnetism_A=bdg_uni._hartree_iterations[0]-bdg_uni._hartree_iterations[1]
    magnetism_B=bdg_uni._hartree_iterations[2]-bdg_uni._hartree_iterations[3]
    axs[0,0].plot(magnetism_A,c='g',marker=markers[0],markersize=s,label=r'$\langle\hat{M}_A\rangle$', linestyle=linestyle[0])
    axs[0,0].plot(magnetism_B,c='b',marker=markers[1],markersize=s,label=r'$\langle\hat{M}_B\rangle$', linestyle=linestyle[1])

    axs[1,0].plot(bdg_uni._fock_iterations[0],c='r',marker=markers[0],markersize=s,label=r'$\Phi^{AB}_{\uparrow\uparrow}$', linestyle=linestyle[0])
    axs[1,0].plot(bdg_uni._fock_iterations[1],c='g',marker=markers[1],markersize=s,label=r'$\Phi^{AB}_{\downarrow\downarrow}$', linestyle=linestyle[1])
    axs[1,0].plot(bdg_uni._fock_iterations[2],c='b',marker=markers[2],markersize=s,label=r'$\Phi^{AA}_{\uparrow\downarrow}$', linestyle=linestyle[2])
    axs[1,0].plot(bdg_uni._fock_iterations[3],c='k',marker=markers[3],markersize=s,label=r'$\Phi^{BB}_{\uparrow\downarrow}$', linestyle=linestyle[3])

    axs[2,0].plot(bdg_uni._gorkov_iterations[0],c='r',marker=markers[0],markersize=s,label=r'$\Delta^{AB}_{\uparrow\uparrow}$', linestyle=linestyle[0])
    axs[2,0].plot(bdg_uni._gorkov_iterations[1],c='g',marker=markers[1],markersize=s,label=r'$\Delta^{AB}_{\downarrow\downarrow}$', linestyle=linestyle[1])
    axs[2,0].plot(bdg_uni._gorkov_iterations[2],c='b',marker=markers[2],markersize=s,label=r'$\Delta^{AA}_{\uparrow\downarrow}$', linestyle=linestyle[2])
    axs[2,0].plot(bdg_uni._gorkov_iterations[3],c='k',marker=markers[3],markersize=s,label=r'$\Delta^{BB}_{\uparrow\downarrow}$', linestyle=linestyle[3])

    axs[3,0].plot(bdg_uni.free_energy,c='r',marker=markers[4],markersize=s,label=f'Free energy')

    # ax4 = ax.twinx()  # instantiate a second axes that shares the same x-axis
    # color = 'tab:red'
    # ax4.set_ylabel('Free energy', color=color)  # we already handled the x-label with ax1
    # ax4.tick_params(axis='y', labelcolor=color)
    # ax4.plot(bdg_uni.free_energy,c='r',marker=markers[4],markersize=s,label=f'Free energy')

    hartree_A_up=bdg_non.hartree(atom='A',spin='up')[0,1]
    hartree_A_dn=bdg_non.hartree(atom='A',spin='dn')[0,1]
    hartree_B_up=bdg_non.hartree(atom='B',spin='up')[0,1]
    hartree_B_dn=bdg_non.hartree(atom='B',spin='dn')[0,1]
    fock_A_B_up_up=bdg_non.fock(atom_i='A', atom_f='B', spin_i='up',spin_f='up', hop_vector=[0,1])[0,0]
    fock_A_B_dn_dn=bdg_non.fock(atom_i='A', atom_f='B', spin_i='dn',spin_f='dn', hop_vector=[0,1])[0,0]
    fock_A_A_up_dn=bdg_non.fock(atom_i='A', atom_f='A', spin_i='up',spin_f='dn', hop_vector=[0,1])[0,0]
    fock_B_B_up_dn=bdg_non.fock(atom_i='B', atom_f='B', spin_i='up',spin_f='dn', hop_vector=[0,1])[0,0]
    gorkov_A_B_up_up=bdg_non.gorkov(atom_i='A', atom_f='B', spin_i='up',spin_f='up', hop_vector=[0,1])[0,0]
    gorkov_A_B_dn_dn=bdg_non.gorkov(atom_i='A', atom_f='B', spin_i='dn',spin_f='dn', hop_vector=[0,1])[0,0]
    gorkov_A_A_up_dn=bdg_non.gorkov(atom_i='A', atom_f='A', spin_i='up',spin_f='dn', hop_vector=[0,1])[0,0]
    gorkov_B_B_up_dn=bdg_non.gorkov(atom_i='B', atom_f='B', spin_i='up',spin_f='dn', hop_vector=[0,1])[0,0]
    free_energy=bdg_non.free_energy

    linestyle=['solid','dashed','dotted','dashdot']

    # color = 'tab:black'
    axs[0,1].tick_params(axis='y')
    magnetism_A=bdg_non._hartree_iterations[0]-bdg_non._hartree_iterations[1]
    magnetism_B=bdg_non._hartree_iterations[2]-bdg_non._hartree_iterations[3]
    axs[0,1].plot(magnetism_A,c='g',marker=markers[0],markersize=s,label=r'$\langle\hat{M}_A\rangle$', linestyle=linestyle[0])
    axs[0,1].plot(magnetism_B,c='b',marker=markers[1],markersize=s,label=r'$\langle\hat{M}_B\rangle$', linestyle=linestyle[1])

    axs[1,1].plot(bdg_non._fock_iterations[0],c='r',marker=markers[0],markersize=s,label=r'$\Phi^{AB}_{\uparrow\uparrow}$', linestyle=linestyle[0])
    axs[1,1].plot(bdg_non._fock_iterations[1],c='g',marker=markers[1],markersize=s,label=r'$\Phi^{AB}_{\downarrow\downarrow}$', linestyle=linestyle[1])
    axs[1,1].plot(bdg_non._fock_iterations[2],c='b',marker=markers[2],markersize=s,label=r'$\Phi^{AA}_{\uparrow\downarrow}$', linestyle=linestyle[2])
    axs[1,1].plot(bdg_non._fock_iterations[3],c='k',marker=markers[3],markersize=s,label=r'$\Phi^{BB}_{\uparrow\downarrow}$', linestyle=linestyle[3])

    axs[2,1].plot(bdg_non._gorkov_iterations[0],c='r',marker=markers[0],markersize=s,label=r'$\Delta^{AB}_{\uparrow\uparrow}$', linestyle=linestyle[0])
    axs[2,1].plot(bdg_non._gorkov_iterations[1],c='g',marker=markers[1],markersize=s,label=r'$\Delta^{AB}_{\downarrow\downarrow}$', linestyle=linestyle[1])
    axs[2,1].plot(bdg_non._gorkov_iterations[2],c='b',marker=markers[2],markersize=s,label=r'$\Delta^{AA}_{\uparrow\downarrow}$', linestyle=linestyle[2])
    axs[2,1].plot(bdg_non._gorkov_iterations[3],c='k',marker=markers[3],markersize=s,label=r'$\Delta^{BB}_{\uparrow\downarrow}$', linestyle=linestyle[3])

    axs[3,1].plot(bdg_non.free_energy,c='r',marker=markers[4],markersize=s,label=f'Free energy')



    # plt.plot(np.real(bdg_uni._gorkov_iterations[1]))
    fig.set_size_inches(w=LATEX_WIDTH, h=1.2*LATEX_WIDTH) 
    title=r'2D multiorbital superconductivity'    

    fig.suptitle(title,x=0.6)
    axs[0,0].set_title(r'Non-unitary triplet')
    axs[0,1].set_title(r'Unitary triplet')
    fig.supxlabel(r'Iterations')
    axs[0,0].set_ylabel(r'Magnetism')
    axs[1,0].set_ylabel(r'Fock')
    axs[2,0].set_ylabel(r'Gorkov')
    axs[3,0].set_ylabel(r'Free energy')
    # plt.subplots_adjust(hspace = 1.5)
    
    plt.tight_layout()

    axs[0,1].legend(loc='upper left', bbox_to_anchor=(1.05,1.05))
    axs[0,1].get_legend().set_title(r'Magnetism')
    axs[1,1].legend(loc='upper left', bbox_to_anchor=(1.05,1.4))
    axs[1,1].get_legend().set_title(r'Fock')
    axs[2,1].legend(loc='upper left', bbox_to_anchor=(1.05,1.15))
    axs[2,1].get_legend().set_title(r'Gorkov')



    # place a text box in lower right axs
    props = dict(boxstyle='round', facecolor='w', alpha=0.2)
    textstr=rf'$\mu={mu}$'+r'''
'''+rf'$t_d={td}$'+r'''
'''+rf'$v={v}$'+r'''
'''+rf'$w={w}$'+'''
'''+rf'$U_v={Uv}$'+r'''
'''+rf'$U_s={Us}$'

    axs[-1,-1].text(1.15, 0.85, textstr, transform=axs[-1,-1].transAxes,
        verticalalignment='top', bbox=props)

    FIGNAME=r'2D multiorbital spin-triplet'
    OUTPUT=os.path.join(FIG,FIGNAME)

    plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight")

    with open(OUTPUT+'.txt', 'w') as f:
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
    
    OUTPUT=os.path.join(FIG,FIGNAME)
    if GIF:
        if not os.path.exists(OUTPUT):
            os.makedirs(OUTPUT)

        files=glob.glob(os.path.join(OUTPUT,'*'))
        
        # OUTPUT=os.path.join(OUTPUT,str(len(files)))
        # plt.savefig(OUTPUT+'.png', bbox_inches = "tight")
        return fig, axs

    else:
        plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight")
        with open(OUTPUT+'.txt', 'w') as f:
            f.write(rf'''Renormalisation of the Hartree, Fock and Gorkov fields on a 
        2-dimensional spinless Weyl-SSH lattice with 
attractive intracell Hubbard U, pairing fermions on sites A and B. 
The fields are calculated self-consistently, with absolute convergence factor 
${absolute_convergence_factor}$
and an initial friction ${init_friction}$ and 
friction ${iter_friction}$ over subsequent iterations.
The lattice is ${n_cells}x{n_cells}$ cells squared.
''')

    plt.close()

def plot_phase_diagram_Delta(phase_diagram_Delta,FIGNAME):
    fig, axs = plt.subplots(5,1,sharex=True)
    axs[0] = phase_diagram_Delta.plot_phase_diagram(axs[0],field_index=1)
    axs[1] = phase_diagram_Delta.plot_phase_diagram(axs[1],field_index=2)
    axs[2] = phase_diagram_Delta.plot_phase_diagram(axs[2],field_index=3,absolute=True)
    axs[3] = phase_diagram_Delta.plot_phase_diagram(axs[3],field_index=4,absolute=True)
    axs[-1] = phase_diagram_Delta.plot_phase_diagram(axs[-1],field_index=0)

    # axs[0].set_ylabel(r'$|\langle\hat{M}\rangle|$')
    axs[0].set_ylabel(r'$\phi_A$')
    axs[1].set_ylabel(r'$\phi_B$')
    axs[2].set_ylabel(r'$|\Phi_{AB}|$')
    axs[3].set_ylabel(r'$|\Delta_{AB}|$')
    axs[-1].set_ylabel('Free energy')
    axs[-1].set_xlabel(r'$\mu$')

    handles, labels = axs[0].get_legend_handles_labels()
    handle_list, label_list = [], []
    for handle, label in zip(handles, labels):
        if label not in label_list:
            handle_list.append(handle)
            label_list.append(label)
    axs[0].legend(handle_list, label_list,loc='upper left', bbox_to_anchor=(1.05,1))
    axs[0].get_legend().set_title(r'$U_v$')

    # plt.subplots_adjust(hspace=1)

    fig.set_size_inches(w=LATEX_WIDTH, h=1.2*LATEX_WIDTH) 

    OUTPUT=os.path.join(FIG,FIGNAME)
    if GIF:
        if not os.path.exists(OUTPUT):
            os.makedirs(OUTPUT)

        files=glob.glob(os.path.join(OUTPUT,'*'))
        
        OUTPUT=os.path.join(OUTPUT,str(len(files)))
        plt.savefig(OUTPUT+'.png', bbox_inches = "tight")

    else:
        plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight")
        with open(OUTPUT+'.txt', 'w') as f:
            f.write(rf'''Self-consistent Hartree, Fock and Gorkov fields on an 
${n_cells}x{n_cells}$ lattice as a function of Coulomb repulsion U for 
various chemical potentials $\mu$.
The absolute convergence factor is
${absolute_convergence_factor}$
and an initial friction ${init_friction}$ and 
friction ${iter_friction}$ over subsequent iterations.
''')
    plt.close()


#############################################################################
################################# Main ######################################
#############################################################################
mu=0.23
s=0
v=0.6
w=1.2
td=0.9
Uv=3.3
Us=5.6
rho=-0.6
rho_shift=2.5
phi=0
chi=1.828633
chi_shift=0.5
n_cells=5

friction=0.95
absolute_convergence_factor=0.00001

bulk_calculation=True

_print=False
bdg = model(mu,Us,Uv,s)
# bdg._set_tightbinding_ham()
# bdg._set_mean_field_hamiltonian()
# bdg._set_hubbard_indices()
# print(bdg._fock)
# bdg.self_consistent_calculation(friction=friction, max_iterations=200, absolute_convergence_factor=absolute_convergence_factor)
# plot_iterations(bdg)

def main():
    bdg_uni = model(mu,Us,Uv,s)
    bdg_uni.self_consistent_calculation(friction=friction, max_iterations=200, absolute_convergence_factor=absolute_convergence_factor)
    FIGNAME=r'2D multiorbital spin-triplet theory'

    chi_shift=0
    bdg_non = model(mu,Us,Uv,s)
    bdg_non.self_consistent_calculation(friction=friction, max_iterations=200, absolute_convergence_factor=absolute_convergence_factor)
    FIGNAME=r'2D multiorbital spin-triplet theory'
    with open(DATA+'unitary.npz', 'wb') as f:
        cPickle.dump(bdg_uni, f)
    with open(DATA+'non_unitary.npz', 'wb') as f:
        cPickle.dump(bdg_non, f)
# main()
bdg_uni=np.load(DATA+'unitary.npz', allow_pickle=True)
bdg_non=np.load(DATA+'non_unitary.npz', allow_pickle=True)

plot_iterations_uni_vs_non(bdg_uni,bdg_non)
exit()

# Phase diagram Uv, mu
rho_shift=2.443
rho=-1.53378
phi=0 #0.6427
chi=0 #1.828633
init_friction=0
iter_friction=0
init_max_iterations=200
iter_max_iterations=200
absolute_convergence_factor=0.00001

GIF=True
GIF=False
plot_initial_renormalisation=None

_print=False
Uvs=[-4.2,-2.8,-1.4,0,1.4,2.8,4.2]
muu=np.arange(-6.1,10.1,0.5)[::-1]
muuA=muu[:int(len(muu)/2)]
muuB=muu[int(len(muu)/2):]
muuu=[muuA,muuB]
data=['A','B']
for i in range(2):
    muu=muuu[i]
    datalabel=data[i]
    rho=-2.4
    rho_shift=0#.443
    phi=0
    chi=0.8
    phase_diagram_Delta=PhaseDiagram(model_Uv)
    phase_diagram_Delta.directory=DATA_Delta
    phase_diagram_Delta.filename=filename
    phase_diagram_Delta.initial_name='initial_convergence_Delta'
    phase_diagram_Delta.initial_title=f'Stoner theory with nonzero Gorkov and Fock pairing\n$\mu={muu[0]:.2f},\, U={Uvs[0]:.2f},\, \,$'
    FIGNAME='Self-consistent_convergence_Delta'
    # phase_diagram_Delta.phase_diagram(muu,dependent_variables,Uvs,init_friction=init_friction,iter_friction=iter_friction,init_max_iterations=init_max_iterations,iter_max_iterations=iter_max_iterations,absolute_convergence_factor=absolute_convergence_factor,initial_renormalisation_plot_function=plot_initial_renormalisation,plot_phase_diagram_function=plot_phase_diagram_Delta,FIGNAME=FIGNAME,datalabel=datalabel)
    plot_phase_diagram_Delta(phase_diagram_Delta,FIGNAME=FIGNAME)
