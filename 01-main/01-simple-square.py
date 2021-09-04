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
################# Simple square model ###################
#########################################################

if n_orbs==1:   
    orbitals=[[0,0,0]]
if n_orbs==2:
    orbitals=[[0,0,0],[0.35,0,0]]

if n_orbs==1:
    hoppings=[[t, 0, 0, [1,0,0]], [t, 0, 0, [0,1,0]]]
if n_orbs==2:
    hoppings=[[t, 0, 1, [0,0,0]],[t, 1, 0, [1,0,0]],[t, 1, 1, [0,1,0]],[t, 0, 0, [0,1,0]]]

basis=[[1,0,0],[0,1,0]]
pbc=[True,True]
#pbc=[False,False]
SO_spin=np.eye(n_spins)
SO_orbital=np.eye(n_orbs)
SO_tensor=-mu*kron(SO_spin,SO_orbital)
############ hopping ##############
hoppings = [[t, 0, 0, [1,0]],
            [t, 0, 0, [0,1]]]
if len(dimensions)==3:
    pbc=[True,True,True]
    #pbc=[False,False,False]
    if dimensions[2]>1:
        basis=[[1,0,0],[0,1,0],[0,0,1]]
        hoppings = [[t, 0, 0, [1,0,0]],
                    [t, 0, 0, [0,1,0]],
                    [t, 0, 0, [0,0,1]]]
############# impurities ##############
impurity_loc = [[0,0,0]]
impurity_spin = np.eye(n_spins)
impurity_orb = np.eye(n_orbs)
impurities = [V, impurity_loc, impurity_spin, impurity_orb]
############# interaction #############
if model=='sc':
    if hubbard=='onsite':
        U_orbs = np.eye(n_orbs)
        U_spin = Pauli_x
        hubbard_SO = -U*kron(U_spin, U_orbs)
        hubbard_U =[[hubbard_SO, [0,0,0]]]
    if hubbard=='nn':
        U_orbs = np.eye(n_orbs)
        U_spin = Pauli_x
        hubbard_SO = -U*kron(U_spin, U_orbs)
        U_orbs_nn = np.eye(n_orbs)
        U_spin_nn = np.eye(n_spins)#+Pauli_x
        hubbard_SO_nn = -U*kron(U_spin_nn, U_orbs_nn)
        hubbard_U =[[hubbard_SO, [0,0,0]],
                    [hubbard_SO_nn, [1,0,0]],
                    [hubbard_SO_nn, [0,1,0]]]
        if len(dimensions)==3:
            if dimensions[2]>1:
                hubbard_U =[[hubbard_SO, [0,0,0]],
                            [hubbard_SO_nn, [1,0,0]],
                            [hubbard_SO_nn, [0,1,0]],
                            [hubbard_SO_nn, [0,0,1]]]
    ############ Mean-field #############
    phiUp, phiDown, Delta= phi_initial, phi_initial, Delta_initial
    varsigma = 0#.035
    Delta1, Delta2 = np.array([varsigma, -varsigma]) + Delta

    hartree_SO = kron(np.diag([phiUp,phiDown]),np.eye(n_orbs))
    initial_hartree = kron(hartree_SO, np.eye(np.prod(dimensions)))
    initial_hartree = np.diag(initial_hartree)

    initial_fock = np.kron(np.kron(Pauli_x,np.eye(n_orbs)),np.ones([np.prod(dimensions),np.prod(dimensions)]))

    gorkov_SO = kron(np.array([[0,+Delta1],[-Delta2,0]]), np.eye(n_orbs))
    initial_gorkov = kron(gorkov_SO, np.eye(np.prod(dimensions)))
############## Energy ###############
epsilon=0.1
increment=0.5
max_val=0 #10*t 
omegas = np.arange(start=-max_val, stop=max_val+increment, step=increment, dtype=complex)+1.0j*epsilon
##########################################################
######################### Main ###########################
##########################################################
if model=='tb':
    model = TB(dimensions, n_spins, basis, SO_tensor, orbitals, hoppings, impurities, epsilon, omegas, pbc)
    model.solve()

elif model=='mf':
    model = BdG(dimensions, n_spins, basis, SO_tensor, orbitals, hoppings, impurities, epsilon, omegas, hartree, fock, gorkov)
    model.solve()

elif model=='sc':
    model = BdG_SC(dimensions, n_spins, basis, SO_tensor, orbitals, hoppings, impurities, epsilon, omegas, initial_hartree, initial_fock, initial_gorkov, hubbard_U, external_hartree=None, external_fock=None, external_gorkov=None,  T=T, friction=friction, max_iterations=max_iterations, eps=eps)

    index1=model.index[0,0,0,0,0]
    index2=model.index[0,0,0,1,0]
    index3=model.index[1,0,0,0,0]

    if hubbard=='nn':
        hartree_indices = [index1, index2]
        fock_indices = [[index1, index3], [index3, index1]]
        gorkov_indices = [[index1,index2], [index2,index1], [index1,index3], [index3,index1]]
    if hubbard=='onsite':
        hartree_indices = [index1, index2]
        fock_indices = [[index1,index2],[index2,index1]]
        gorkov_indices = [[index1,index2], [index2,index1]]

    model.set_field_tracking(hartree_indices=hartree_indices, fock_indices=fock_indices, gorkov_indices=gorkov_indices)

    hartree, fock, gorkov, iterations = model.solve(dos=True)
    
#    y1=model.hartree_list[0]
#    y2=model.fock_list[0]
#    y3=model.gorkov_list[0]
#    y4=model.gorkov_list[2]
#    y5=model.gorkov_list[4]
    
    x0,y0,z0=0,0,0
    x1,y1,z1=1,0,0
    s,t=0,1
    m,n=0,0
    
    hartree=model.hartree()
    print('Hartree onsite:')
    print(hartree[:,y0,z0,s,m])

    fock=model.fock()
    renormalised_t=[fock[(x0+i)%n_x,y0,z0,s,m,(x1+i)%n_x,y0,z0,s,m] for i in range(n_x)]
    print('Renormalised t:')
    print(renormalised_t)

    gorkov=model.gorkov()
    g_onsite=[gorkov[(x0+i)%n_x,y0,z0,s,m,(x0+i)%n_x,y0,z0,t,n] for i in range(n_x)]
    g_nn=[gorkov[(x0+i)%n_x,y0,z0,s,m,(x1+i)%n_x,y0,z0,s,m] for i in range(n_x)]
    print('Gorkov onsite:')
    print(g_onsite)
    print('Gorkov nn:')
    print(g_nn)


    index  = [index1,index2]
    ind = np.argwhere(((model.hubbard_indices)[0]==index[0], (model.hubbard_indices)[1]==index[1]))

    index = [index1,index2]
    ind = np.argwhere(np.logical_and((model.hubbard_indices)[0]==index[0], (model.hubbard_indices)[1]==index[1]))

    if hubbard=='nn':

        index  = [index1,index3]
        ind = np.argwhere(np.logical_and((model.hubbard_indices)[0]==index[0], (model.hubbard_indices)[1]==index[1]))

        index = [index1,index3]
        ind = np.argwhere(np.logical_and((model.hubbard_indices)[0]==index[0], (model.hubbard_indices)[1]==index[1]))

    print('Free energy')
    print(model.free_energy)


#     import matplotlib.pyplot as plt
#     from matplotlib.ticker import MaxNLocator
# 
#     ax = plt.figure().gca()
#     plt.plot(y1,color='orange',marker='+')
#     plt.plot(y2,color='green',marker='x')
#     plt.plot(y3,color='blue',marker='.')
#     plt.plot(y4,color='red',marker='o')
#     ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#     plt.show()
#     plt.close()

with open(confname+'.npz', 'wb') as f:
    cPickle.dump(model, f)
########################## temp ###########################
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)
#fig.set_size_inches(w=latex_width, h=4.5) 

omega=0
layer=(slice(None),slice(None),0)
ldos = LDOS(model.density_of_states, omegas, omega, trace_over=True)
ldos = LocalDensityOfStates(fig, ax, model, ldos, omega)
fig=ldos.plot_3D()
fig.show()
fig,ax=ldos.imshow(layer)
plt.show()

#######
exit()
#######
ham = model.set_tb_hamiltonian()
w,v = la.eigh(ham, overwrite_a=True)
dm = model._density_matrix(v) 
model.density_of_states = DOS(model.omegas,w,dm)
unshaped=model.density_of_states
model.density_of_states = model._density_of_states(w,v)
print('\nHamiltonian:')
print('Unshaped:')
ham=np.real(ham)
print(np.diag(ham))
print(np.shape(ham))
dim = np.ravel(np.array([model.extended_dimensions,model.extended_dimensions]).T)
shaped=np.reshape(ham, dim)
print('Reshaped:')
print(shaped[:,:,0,0,0])
print(np.shape(ham))
n_dof=model.n_dof
unshaped=np.reshape(ham, [n_dof, n_dof])
print('Check')
print(np.array_equal(ham, unshaped))
print('------------------')
print('\nDensity of states')
print('Unshaped:')
print(np.shape(model.density_of_states))
dim = np.ravel(np.array([model.extended_dimensions,model.extended_dimensions]).T)
shaped=np.reshape(ham, dim)
print('Reshaped:')
print(np.shape(ham))
n_dof=model.n_dof
unshaped=np.reshape(ham, [n_dof, n_dof])
print('Check')
print(np.array_equal(ham, unshaped))

omegas=self.omegas
dimensions=np.append(self.extended_dimensions,len(omegas))

density_matrix=self._density_matrix(v)

density_of_states=DOS(omegas,w,density_matrix)

self.density_of_states = np.reshape(density_of_states, dimensions)

temp = np.copy(self.density_of_states)

self.density_of_states = np.array(
        [[[[[density_of_states[index[x, y, z, s, m]]
        for m in range(dimensions[4])]
        for s in range(dimensions[3])]
        for z in range(dimensions[2])]
        for y in range(dimensions[1])] 
        for x in range(dimensions[0])])
print(temp[:,0,0,0,0,0])
print(model.density_of_states[:,0,0,0,0,0])
print(np.array_equal(temp,self.density_of_states))
