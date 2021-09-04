from lib import *
#########################################################
################# from .conf import * ###################
#########################################################
if len(sys.argv) < 2:
        print('Please supply conf file')
        sys.exit()
    
conf = sys.argv[-1]

confname = conf.split('.conf')[0]
config_module = import_path(os.path.join(MAIN,conf))
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

confname = os.path.basename(confname) 
fileName = os.path.dirname(sys.argv[1])
fileName = os.path.basename(fileName)

confname = os.path.join(DATA,fileName,confname)
#########################################################
##################### Tightbinding ######################
#########################################################
def tb():
    global dos, exec_time
    tb_model = TB(dimensions, n_spins, basis, SO_tensor, orbitals, hoppings, impurities, epsilon, omegas, silent=True)
    dos = tb_model.solve(silent=True)
    exec_time = tb_model.exec_time
    # mem = tb_model.mem
#########################################################
############ Bogoliubov-de Gennes mean field ############
#########################################################
def bdg():
    bdg_model = BdG()
    dos, ados = tb_model.solve(silent=True)
    exec_time = bdg_model.exec_time
    mem =bdg_model.mem

    with open(confname+'.npz', 'wb') as f:
        np.savez(f,
        dos=dos,
        ados=ados,
        exec_time=exec_time,
        mem=mem)
#########################################################
########## Bogoliubov-de Gennes self-consistent #########
#########################################################
def bdg_singlet_sc():
    bdg = BdG_SC(dimensions, n_spins, basis, SO_tensor, orbitals, hoppings, impurities, epsilon, omegas, initial_hartree, initial_fock, initial_gorkov, hubbard_U, external_hartree=None, external_fock=None, external_gorkov=None,  T=T, friction=friction, max_iterations=max_iterations, eps=eps, silent=False)

    index1=bdg.index[0,0,0,0,0]
    index2=bdg.index[0,0,0,1,0]

    hartree_indices = [index1, index2]
    gorkov_indices = [[index1,index2], [index2,index1]]

    bdg.set_field_tracking(hartree_indices=hartree_indices, gorkov_indices=gorkov_indices)

    iteration = iter(bdg)

    hartree, fock, gorkov, iterations = bdg.self_consistent()

    hartree_list=bdg.hartree_list
    fock_list=[]
    gorkov_list=bdg.gorkov_list

    # dos=bdg.dos
    # ados=bdg.ados

    exec_time = bdg.exec_time
    mem = bdg.mem
    with open(confname+'.npz', 'wb') as f:
        np.savez(f,
        # dos=dos,
        # ados=ados,
        hartree=hartree_list,
        fock=fock_list,
        gorkov=gorkov_list,
        exec_time=exec_time,
        iterations=iterations,
        mem=mem)
#########################################################
######################## Main ###########################
#########################################################
if model=='tb': 
    mem = memory_usage(tb)

    with open(confname+'.npz', 'wb') as f:
        np.savez_compressed(f,
        dos=dos,
        exec_time=exec_time,
        mem=mem)
if model=='bdg': bdg()
if model=='bdg_singlet_sc': bdg_singlet_sc()
