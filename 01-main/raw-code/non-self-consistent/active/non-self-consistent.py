import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
import sys
sys.path.append('../../main/')
from lib import *
if len(sys.argv) < 2:
	print('Please supply name of configuration file')
	sys.exit()
    
out_folder='../out/'

configuration = sys.argv[1]
confname = configuration.split('.conf')[0]
config_module = import_path(configuration)    
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

index=Index(dof, n_x, n_y, n_sites,n_spins,n_orbitals)
########################################################################
############################### Main ###################################
########################################################################
t0=time.perf_counter()
h0 = H0(onsite_tensor, nn_tensor_x, nn_tensor_y, impurity_tensor,impurity_locations, n_x, n_y)
mf = Mean_field(SC_tensor, n_sites)
ham = H(h0, mf)
del(h0)
del(mf)
gc.collect() # Garbage collector prevents MemoryError
memory=sys.getsizeof(ham)
logging.debug('matrix size (MiB):')
logging.debug(memory/1048576)
logging.debug('running la.eigh')
w,v = la.eigh(ham, overwrite_a=True)

logging.debug('running DOS')
density_matrix = DM(v, n_x,n_y,index,n_sites,n_spins,n_orbitals)
dos = DOS(omegas, w, density_matrix)
density_matrix = ADM(v, n_x,n_y,index,n_sites,n_spins,n_orbitals)
ados = DOS(omegas, w, density_matrix)

t1=time.perf_counter()
executionTime=t1-t0
logging.debug('time taken (s):')
logging.debug(executionTime)
tt=time.perf_counter()

logging.debug('saving')
with open(out_folder+confname+'.npz', 'wb') as f:
    np.savez(f,
    dos=dos,
    ados=ados,
    executionTime=executionTime,
    memory=memory)