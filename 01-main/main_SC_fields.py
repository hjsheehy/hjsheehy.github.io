from lib import *
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
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
########################################################
####################### Main ###########################
########################################################
bdg = BdG_SC(dimensions, n_spins, basis, SO_tensor, orbitals, hoppings, impurities, epsilon, omegas, initial_hartree, initial_fock, initial_gorkov, hubbard_U, external_hartree=None, external_fock=None, external_gorkov=None,  T=T, friction=friction, max_iterations=max_iterations, eps=eps, silent=False)

index1=bdg.index[0,0,0,0,0]
index2=bdg.index[0,0,0,1,0]

hartree_indices = [index1, index2]
gorkov_indices = [[index1,index2], [index2,index1]]

bdg.set_field_tracking(hartree_indices=hartree_indices, gorkov_indices=gorkov_indices)

iteration = iter(bdg)

hartree, fock, gorkov, iterations = bdg.self_consistent()

y1=bdg.gorkov_list[0]
y2=-bdg.gorkov_list[1]
y3=bdg.hartree_list[0]
y4=bdg.hartree_list[1]

print(y3)
print(y1)

ax = plt.figure().gca()
plt.plot(y1,color='orange',marker='+')
plt.plot(y2,color='green',marker='x')
plt.plot(y3,color='blue',marker='.')
plt.plot(y4,color='red',marker='o')
# plt.ylim([0,2])
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.show()
# plt.close()
# print(GorkovUp_list[-1])
# print(phiUp_list[-1])
# print(len(GorkovUp_list))
# print(executionTime)

omega=0
layer=0
bdg.density_of_states()
ldos = LDOS(bdg.dos, bdg.omegas, omega, trace_over=True, layer=layer)
x=np.fft.fftshift(ldos)
plt.imshow(x, extent=[-x.shape[1]/2., x.shape[1]/2., -x.shape[0]/2., x.shape[0]/2. ])
plt.show()
