from lib import *
#########################################################
################# from .conf import * ###################
#########################################################
if len(sys.argv) < 2:
        print('Please supply conf file')
        sys.exit()
    
conf = sys.argv[-1]

confName = conf.split('.conf')[0]
config_module = import_path(os.path.join(MAIN,conf))
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

confName = os.path.basename(confName) 
fileName = os.path.dirname(sys.argv[1])
fileName = os.path.basename(fileName)
confName = os.path.join(DATA,fileName,confName)

#########################################################
################# Simple square model ###################
#########################################################
temperature=0
max_iterations=10000
friction=0.9
absolute_convergence_factor=0.00001
############## Energy ###############
resolution=0.1
increment=0.01
max_val=10
energy_interval = np.arange(start=-max_val, stop=max_val+increment, step=increment)

def Main(mu, Delta_v):
    ##########################################################
    ######################### Main ###########################
    ##########################################################
    model = bogoliubov_de_gennes(dimensions=[n,1], n_spins=2, basis=[[1,0],[0,1]], orbitals=[[0.35,0.5]], pbc=[True,True])
    model.set_onsite(-mu)
    model.set_hopping(-t, 0, 0, [1,0],'$-t$')
    model.set_hubbard_u(-U*Pauli_x, orb_i=0, orb_f=0, hop_vector=[0,0])

    # model = bogoliubov_de_gennes(dimensions=[n,1], n_spins=2, basis=[[1,0],[0,1]], orbitals=[[0.35,0.5],[0.65,0.5]], pbc=[True,True])
    # model.set_onsite(-mu)
    # model.set_hopping(-v, 0, 1, [0,0],'$-v$')
    # model.set_hopping(-w, 1, 0, [1,0],'$-w$')
    # model.set_hopping(-w, 0, 0, [0,1],'$-w$')
    # model.set_hopping(-w, 1, 1, [0,1],'$-w$')
    # model.set_impurities(V, locations=[0],spins=[0,1],orbitals=[0,1])

    # model.set_hubbard_u(-U_T*np.eye(2), orb_i=0, orb_f=1, hop_vector=[0,0])
    # hubbard_R=U_R*np.eye(2)
    # model.set_hubbard_u(hubbard_R, orb_i=1, orb_f=0, hop_vector=[1,0])
    # model.set_hubbard_u(hubbard_R, orb_i=0, orb_f=0, hop_vector=[0,1])
    # model.set_hubbard_u(hubbard_R, orb_i=1, orb_f=1, hop_vector=[0,1])

    spin_sector=np.eye(2)
    spin_sector[[1,1]]=0

    def set_mean_fields():
        
        if state=='BCS':
            model.set_hartree(phi)
            model.set_gorkov(Delta_v*1.0j*Pauli_y, orb_i=0, orb_f=0, hop_vector=[0,0])

        if state=='INT':
            model.set_hartree([1.25*phi,0.75*phi])

            model.set_fock(zeta_v*spin_sector, orb_i=0, orb_f=1, hop_vector=[0,0])
            model.set_fock(zeta_w*np.eye(2), orb_i=1, orb_f=0, hop_vector=[1,0])
            model.set_fock(zeta_w*np.eye(2), orb_i=0, orb_f=0, hop_vector=[0,1])
            model.set_fock(zeta_w*np.eye(2), orb_i=1, orb_f=1, hop_vector=[0,1])

            model.set_gorkov(Delta_v*spin_sector, orb_i=0, orb_f=1, hop_vector=[0,0])
            model.set_gorkov(Delta_w*np.eye(2), orb_i=1, orb_f=0, hop_vector=[1,0])
            model.set_gorkov(Delta_w*np.eye(2), orb_i=0, orb_f=0, hop_vector=[0,1])
            model.set_gorkov(Delta_w*np.eye(2), orb_i=1, orb_f=1, hop_vector=[0,1])

        if state=='unitary':
            model.set_hartree([1.25*phi,0.75*phi])

            model.set_fock(zeta_v*np.eye(2), orb_i=0, orb_f=1, hop_vector=[0,0])
            model.set_fock(zeta_w*np.eye(2), orb_i=1, orb_f=0, hop_vector=[1,0])
            model.set_fock(zeta_w*np.eye(2), orb_i=0, orb_f=0, hop_vector=[0,1])
            model.set_fock(zeta_w*np.eye(2), orb_i=1, orb_f=1, hop_vector=[0,1])

            model.set_gorkov(Delta_v*np.eye(2), orb_i=0, orb_f=1, hop_vector=[0,0])
            model.set_gorkov(Delta_w*np.eye(2), orb_i=1, orb_f=0, hop_vector=[1,0])
            model.set_gorkov(Delta_w*np.eye(2), orb_i=0, orb_f=0, hop_vector=[0,1])
            model.set_gorkov(Delta_w*np.eye(2), orb_i=1, orb_f=1, hop_vector=[0,1])

    set_mean_fields()

    _print=False
    _print=True

    model.record_gorkov(location_a=[0,0], location_b=[0,0], spin_a=0, spin_b=1, orbital_a=0, orbital_b=0, _print=_print)
    model.record_hartree([0,0], 0, 0, _print)
    # model.record_hartree([0,0], 1, 1, _print)
    # model.record_gorkov(location_a=[5,0], location_b=[5,0], spin_a=0, spin_b=0, orbital_a=0, orbital_b=1, _print=_print)
    # model.record_gorkov(location_a=[4,0], location_b=[5,0], spin_a=0, spin_b=0, orbital_a=1, orbital_b=0, _print=_print)

    model.set_max_iterations(max_iterations)
    model.set_friction(friction)
    model.set_absolute_convergence_factor(absolute_convergence_factor)
    model.set_temperature(temperature)
    eigenvalues,eigenvectors=model.self_consistent_calculation()    
            
    # Trash unnecessary data
    del(model._tb_ham)
    del(model._hubbard_u)

    data=processed_data(model, energy_interval, resolution, eigenvalues, eigenvectors)
    print(np.max(data.gorkov()))
    x=y=m=0
    return model.free_energy

x=[]
y=[]
Main(mu,Delta_v)
# exit()
max_iterations=1
friction=1

# for mu in np.linspace(-5,5,41):
#     x.append(mu)
#     y.append(Main(mu, Delta_v))
for delta_v in np.linspace(-2.5,2.5,41):
    x.append(delta_v)
    y.append(Main(mu,delta_v))

plt.plot(x,y)
plt.scatter(x,y)
plt.show()

exit()
########################################################
########################## Plot ########################
########################################################
latex_width=4.7747
n_pts=40    
def fig_a():
    energy=0

    fig, axs = plt.subplots(1)
    plt_data = plot_data(fig,axs,data)

    fig, axs = plt_data.differential_current_map(energy, layer=(ALL,ALL), orbital=[0,1], spin=None, cartesian=True, n_pts=n_pts, gaussian_mean=0.1)

    label_x = r'$x/a_x$'
    label_y = r'$y/a_y$'
    axs.set_ylabel(label_y)
    axs.set_xlabel(label_x)

    plt.tight_layout()

    axs.set_title('zero-bias local density of states')
    return fig

y=data.centre[1]
gorkov=data.gorkov()
Delta_T=np.real([[gorkov[x,y,0,s,x,y,1,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
Delta_T=np.sum(Delta_T,0)/data.n_spins
Delta_R=np.real([[gorkov[x-1,y,1,s,x,y,0,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
Delta_R=np.sum(Delta_R,0)/data.n_spins
print(f'Delta_T={Delta_T[-1]}')
print(f'Delta_R={Delta_R[-1]}')

fig_a()
plt.show()
