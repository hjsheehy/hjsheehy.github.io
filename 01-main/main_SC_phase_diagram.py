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
xx=0
yy=0
zz=0
s=0
t=1
mm=0

# phiUp_list=[bdg.hartree[xx,yy,zz,s,mm]]
# phiDown_list=[bdg.hartree[xx,yy,zz,t,mm]]
# GorkovUp_list=[bdg.gorkov[xx,yy,zz,s,t,mm,mm]]
# GorkovDown_list=[-bdg.gorkov[xx,yy,zz,t,s,mm,mm]]
#########################################################
######################## Main ###########################
#########################################################
#########################################################
################# Plotting -delete me ###################
#########################################################
######################### LDOS ##########################
# ldos = np.einsum('xyzssmme->xyze', dos)

# layer=0
# omega=0
# index = FindNearestValueOfArray(np.real(omegas), omega)
# result = ldos[:,:,layer,index]

# plt.imshow(result)
# plt.show()
####################### Phase diagram ####################
n_U=4
n_T=60
U_list = np.linspace(1.2, 2.8, n_U)
#U_list=[2.8]
T_list = np.linspace(0, 0.5, n_T)

# Fock_data = np.zeros([n_T,n_U,dimensions[0],dimensions[1],dimensions[2],n_U-1,n_spins,n_spins,n_orbs,n_orbs])

# Hartree_data = np.zeros([n_T,n_U,dimensions[0],dimensions[1],dimensions[2],n_spins,n_orbs])

# Gorkov_data= np.zeros([n_T,n_U,dimensions[0],dimensions[1],dimensions[2],n_U-1,n_spins,n_spins,n_orbs,n_orbs])

bdg = BdG_SC(dimensions, n_spins, basis, SO_tensor, orbitals, hoppings, impurities, epsilon, omegas, initial_hartree, initial_fock, initial_gorkov, hubbard_U, external_hartree=None, external_fock=None, external_gorkov=None,  T=T, friction=friction, max_iterations=max_iterations, eps=eps, silent=False)

index1=bdg.index[0,0,0,0,0]
index2=bdg.index[0,0,0,1,0]

Hartree_data=np.zeros([n_T,n_U,bdg.n_dof])
Fock_data=np.zeros([n_T,n_U,bdg.n_dof,bdg.n_dof])
Gorkov_data=np.zeros([n_T,n_U,bdg.n_dof,bdg.n_dof])

initial_hartree_,initial_gorkov_,initial_fock_ = initial_hartree,initial_gorkov,initial_fock

t = time.time()
for j in tqdm(range(n_U)):
    U=U_list[j]
    initial_hartree,initial_gorkov,initial_fock = np.copy(initial_hartree_),np.copy(initial_gorkov_),np.copy(initial_fock_)
    for i in tqdm(range(n_T),leave=False):
        T=T_list[i]
        
        hubbard_SO = -U*np.kron(U_spin, U_orbs)
        hubbard_U=np.kron(hubbard_SO, np.eye(n_x*n_y))

        bdg = BdG_SC(dimensions, n_spins, basis, SO_tensor, orbitals, hoppings, impurities, epsilon, omegas, initial_hartree, initial_fock, initial_gorkov, hubbard_U, external_hartree=None, external_fock=None, external_gorkov=None,  T=T, friction=friction, max_iterations=max_iterations, eps=eps, silent=False)
        hartree, fock, gorkov, iterations = bdg.self_consistent()

        Hartree_data[i,j]=hartree
        Fock_data[i,j]=fock
        Gorkov_data[i,j]=gorkov

        initial_hartree,initial_fock,initial_gorkov = hartree,fock,gorkov
        
        # print(bdg.iterations)
        # i+=1
        
        # percent_complete = 100*(i+j*len(T_list))/(len(U_list)*len(T_list))
        # loading=str('|'+int(percent_complete/10)*'█'+int(10-percent_complete/10)*'_'+'|')
        # print(f"{loading} {percent_complete}%")
    j+=1
exec_time = time.time() - t
print(exec_time)

# with open(out_folder+confname+'.npz', 'wb') as f:
#     np.savez(f,Hartree_data=Hartree_data,Fock_data=Fock_data,Gorkov_data=Gorkov_data,U_list=U_list,T_list=T_list)

x=T_list
xx=0
yy=0
zz=0
s=0
t=1
mm=0
u=0
y=[]
color=['orange','green','blue','red']
marker=['+','x','.','o']
ax = plt.figure().gca()
for i in range(n_U):
    y.append(Gorkov_data[:,i,index1,index2])
    plt.plot(x,y[i],color=color[i],marker=marker[i])

# plt.ylim([0,0.25])
# ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.show()
