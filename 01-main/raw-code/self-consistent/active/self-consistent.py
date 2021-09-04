import sys
sys.path.append('../../main/')
from lib import *

if len(sys.argv) < 2:
	print('Please supply conf file')
	sys.exit()
    
out_folder='../out/'
    
conf = sys.argv[1]
confname = conf.split('.conf')[0]
config_module = import_path(conf)
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

index=Index(dof,n_x,n_y,n_sites,n_spins,n_orbitals)
Hubbard_U = indexed(Hubbard_tensor,n_x,n_y,index,n_sites,n_spins,n_orbitals)
if 'Gorkov_tensor' in globals():
    Gorkov = indexed(Gorkov_tensor, n_x,n_y,index,n_sites,n_spins,n_orbitals)

def Update(Hartree, Fock, Gorkov):  
    global w, v
    h0 = H0(onsite_tensor, nn_tensor_x, nn_tensor_y, impurity_tensor, impurity_locations, n_x, n_y)
    ham = H_BdG(h0, Hartree, Fock, Gorkov, index, n_x,n_y,n_sites,n_spins,n_orbitals,dof)
    del(h0)
    w,v = la.eigh(ham, overwrite_a=True)

    dm=DM(v, n_x,n_y,index,n_sites,n_spins,n_orbitals)
    dm=Normal_Fermi(dm, dof, w, T)
    adm=ADM(v, n_x,n_y,index,n_sites,n_spins,n_orbitals)
    adm=Anomalous_Fermi(adm, dof, w, T)

    Hartree_update, Fock_update, Gorkov_update = HartreeFockGorkov(Hubbard_U, dm, adm)

    Hartree_update = (1-friction)*Hartree_update+friction*Hartree
    Fock_update = (1-friction)*Fock_update+friction*Fock
    Gorkov_update = (1-friction)*Gorkov_update+friction*Gorkov
    return Hartree_update, Fock_update, Gorkov_update
##########################################################
########################## Main ##########################
########################################################## 
def Main(Hartree, Fock, Gorkov):
    for i in range(max_iterations):
        gc.collect() # Garbage collector prevents MemError
        Hartree_update, Fock_update, Gorkov_update = Update(Hartree, Fock, Gorkov) #remove 0 later
        if (np.allclose(Hartree,Hartree_update,atol=eps) & np.allclose(Fock,Fock_update,atol=eps) & np.allclose(Gorkov,Gorkov_update,atol=eps)):
            return Hartree_update, Fock_update, Gorkov_update, i+1
        Hartree, Fock, Gorkov = Hartree_update, Fock_update, Gorkov_update
    Hartree_update, Fock_update, Gorkov_update = float('nan')*Hartree, float('nan')*Fock, float('nan')*Gorkov
    return Hartree_update, Fock_update, Gorkov_update, max_iterations
 
startTime = time.time()
Hartree, Fock, Gorkov, iterations = Main(Hartree, Fock, Gorkov)
executionTime = (time.time() - startTime)
if Include_DOS==False:
    with open(out_folder+confname+'.npz', 'wb') as f:
        np.savez(f,Hartree=Hartree,Fock=Fock,Gorkov=Gorkov,iterations=iterations,executionTime=executionTime)
        print(iterations)
else:
    print('Running DOS')
    density_matrix = DM(v, n_x,n_y,index,n_sites,n_spins,n_orbitals)
    dos = DOS(omegas, w, density_matrix)
    density_matrix = ADM(v, n_x,n_y,index,n_sites,n_spins,n_orbitals)
    ados = DOS(omegas, w, density_matrix)
    with open(out_folder+confname+'.npz', 'wb') as f:
        np.savez(f,
        Hartree=Hartree,Fock=Fock,Gorkov=Gorkov,
        iterations=iterations,
        executionTime=executionTime,
        dos=dos,
        ados=ados)