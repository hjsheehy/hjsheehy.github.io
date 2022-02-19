from conf import *
#########################################################
################# Simple square bdg ###################
#########################################################
A=Atom([0,0],'A')
B=Atom([0.25,1],'B')
A.add_orbital('s')
B.add_orbital('s')
lattice_vectors=[[1,0],[0,1]]
bdg=BogoliubovdeGennes(lattice_vectors,'SSH')
bdg.add_atom(A)
bdg.add_atom(B)
bdg.n_spins=2
mu=-0.57
w=0.79
v=1
phi=5.32
zeta_v=3.52
zeta_w=2.23
Delta_v=3.92
Delta_w=2.53

n_cells=43

bdg.cut_piece(n_cells, [0])
bdg.set_onsite(-mu,orbital='s')
bdg.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B', label='v')
bdg.set_hopping(-w,hop_vector=[1,0],atom_i='B',atom_f='A', label='w')
bdg.add_impurities(V,[0,0],label='V')
#########################################################
temperature=0
absolute_convergence_factor=0.00001
friction=0.7
max_iterations=1000

bdg.set_temperature(temperature)

bdg.set_hubbard_u(-U_v*np.eye(2), atom_i='A', atom_f='B', hop_vector=[0,0])
hubbard_R=U_w*np.eye(2)
bdg.set_hubbard_u(hubbard_R, atom_i='B', atom_f='A', hop_vector=[1,0])

spin_sector=np.eye(2)
spin_sector[[1,1]]=0

def set_mean_fields(state):

    if state=='INT':
        bdg.set_hartree([1.25*phi,0.75*phi])

        bdg.set_fock(zeta_v*spin_sector, atom_i='A', atom_f='B', hop_vector=[0,0])
        bdg.set_fock(zeta_w*np.eye(2), atom_i='B', atom_f='A', hop_vector=[1,0])

        bdg.set_gorkov(Delta_v*spin_sector, atom_i='A', atom_f='B', hop_vector=[0,0])
        bdg.set_gorkov(Delta_w*np.eye(2), atom_i='B', atom_f='A', hop_vector=[1,0])

    if state=='unitary':
        bdg.set_hartree([1.25*phi,0.75*phi])

        bdg.set_fock(zeta_v*np.eye(2), atom_i='A', atom_f='B', hop_vector=[0,0])
        bdg.set_fock(zeta_w*np.eye(2), atom_i='B', atom_f='A', hop_vector=[1,0])

        bdg.set_gorkov(Delta_v*np.eye(2), atom_i='A', atom_f='B', hop_vector=[0,0])
        bdg.set_gorkov(Delta_w*np.eye(2), atom_i='B', atom_f='A', hop_vector=[1,0])

if config_file:
    
    if layer_no(CONFNAME)==0 and parent_sim==SIM_NAME:
        
        set_mean_fields(state=INITIAL_FIELD)
    else:
        # i.e. if germ or if perturbation from parent_sim

        parent_data = os.path.join(DATA,parent_sim,parent_filename+'.npz')
        parent_data = glob.glob(parent_data)

        # minimum free energy:
        free_energy = [file.split('_')[-1] for file in parent_data]
        min_val = min(free_energy)
        min_index = free_energy.index(min_val)
        parent_filename = parent_data[min_index]
        
        with open(parent_filename, 'rb') as f:
            data = cPickle.load(f)

        bdg._hartree=data._hartree
        bdg._fock=data._fock
        bdg._gorkov=data._gorkov
        bdg._hubbard_indices=data._hubbard_indices
        bdg._anomalous_indices=data._anomalous_indices
        bdg.U_entries=data.U_entries

        # if np.any(np.isnan(hartree)):
        #     print('is nan')
            
        #     exit()

        #     model.reset_hartree()
        #     model.reset_fock()
        #     model.reset_gorkov()

        #     phi=np.random.uniform(low=0.1, high=1.3)
        #     zeta_v=np.random.uniform(low=0.1, high=2.3)
        #     zeta_w=np.random.uniform(low=0.1, high=2.3)
        #     Delta_v=np.random.uniform(low=1.1, high=4.3)
        #     Delta_w=np.random.uniform(low=1.1, high=4.3)
            
        #     set_mean_fields()

bdg.self_consistent_calculation(friction=friction, max_iterations=max_iterations, absolute_convergence_factor=absolute_convergence_factor)

energy_interval=np.linspace(-4,4,51)
resolution=0.1
bdg.calculate_greens_function(energy_interval,resolution)

save_data(bdg_model=bdg, CONFNAME=CONFNAME, SIM_NAME=SIM_NAME)
