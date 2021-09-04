import matplotlib.pyplot as plt #remove later 
import sys
sys.path.append('.')
from lib import *

if len(sys.argv) < 2:
	print('Please supply conf file')
	sys.exit()
    
out_folder='.'
    
conf = sys.argv[1]
confname = conf.split('.conf')[0]
config_module = import_path(conf)
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

index=Index(dimensions, n_sites, n_spins, n_orbitals)

def Main():  
    global w, v
    h0 = H0(lattice_basis, site_vectors, site_tensors,
            hopping_vectors, hopping_tensors, PBC,
            impurity_locations, impurity_tensors,
            n_spins, n_orbitals)
    ham = H(h0, mean_field)
    del(h0)
    w,v = la.eigh(ham, overwrite_a=True)
    dm=DM(v, n, n_spins, n_orbitals, index)
    dos = DOS(omegas, w, dm)
    return dos
##########################################################
########################## Main ##########################
##########################################################
# dos = Main()
# layer=1
# omega=0
# omega_index=FindNearestValueOfArray(omegas,omega)
# result=(np.einsum('xyssooe->xye',dos[:,:,layer]))
# result=result[:,:,omega_index]
# plt.imshow(result)
# plt.show()

print(n_dim)
