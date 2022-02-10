from lib import *
#########################################################
#################### BCS wire model #####################
#########################################################
dimensions=[1,1]
mu=-.7
t=1
n_spins=2
orbitals=[[0,0]]
basis=[[1,0],[0,1]]
pbc=[True]
############# interaction #############
hubbard_U = 2.8*Pauli_x
# hubbard_U = 2.8j*Pauli_y
############ Mean-field #############
phiUp, phiDown, Delta = 1.4,3.4,1.00
############## Energy ###############
resolution=0.1
increment=0.01
max_val=10*t 
energy_interval = np.arange(start=-max_val, stop=max_val+increment, step=increment)
##########################################################
######################### Main ###########################
##########################################################
model = bogoliubov_de_gennes(dimensions, n_spins, basis, orbitals)
model.set_onsite(-mu)
model.set_hopping(-t, 0, 0, [1,0])
model.set_hopping(-t, 0, 0, [0,1])
model.set_impurities(0, [0,0])
model.set_hubbard_u(hubbard_U, None, None, [0])
model.set_hartree([phiUp,phiDown])
model.set_gorkov(Delta*(1.0j*Pauli_y))
_print=False

model.record_hartree([0,0], 0, 0, _print)
model.record_gorkov([0,0], [0,0], 0, 1, 0, 0, _print)

model.set_max_iterations(max_iterations=10000)
model.set_friction(0)
model.set_temperature(0.05)
model.set_absolute_convergence_factor(0.00001)
##########################################################
#################### Genetic algorithm ###################
##########################################################
from geneticalgorithm import geneticalgorithm as ga

model._set_hubbard_indices()

l0=np.shape(model._hubbard_u)[0]
l1=len(model._hubbard_indices[0])

# Random uniform distribution of values on the interval (0, height)
height = 5
height_hartree=height_fock=height_gorkov=height

hartree = height_hartree*np.random.rand(l0)
fock = height_fock*np.random.rand(l1)
gorkov = height_gorkov*np.random.rand(l1)
solution = np.concatenate((hartree, fock, gorkov), axis=0)
l=len(solution)

def fitness_function(solution):

    model.hartree, model.fock, model.gorkov = np.split(solution, [l0,l0+l0])

    ############ Solve mean-field Hamiltonian #############

    model.set_mean_field_hamiltonian()
    w,v = model.solve()

    trace_density, density, anomalous_density = model._thermal_density_matrix(w,v)

    hartree, fock, gorkov = model._set_fields(trace_density, density, anomalous_density)

    ##################### Free energy ######################

    free_energy = model.calculate_free_energy(w, trace_density, density, anomalous_density)

    ####################### Friction ########################

    hartree = (1-model.friction)*hartree + model.friction*model.hartree
    fock = (1-model.friction)*fock + model.friction*model.fock
    gorkov = (1-model.friction)*gorkov + model.friction*model.gorkov

    ###################### Update fields #####################

    model.hartree, model.fock, model.gorkov = hartree, fock, gorkov 

    return free_energy


####################### Parameters #######################
varbound=np.array([[-5,5]]*l)

#################### Initialise ####################
gmodel=ga(function=fitness_function,dimension=l,variable_type='real',variable_boundaries=varbound)

####################### Main ########################
gmodel.run()

################ Collect final data #################
hartree, fock, gorkov = np.split(solution, [l0,l0+l0])
model.hartree, model.fock, model.gorkov = hartree, fock, gorkov
w,v = model.solve()
eigenvalues, eigenvectors = w,v
data=processed_data(model, energy_interval, resolution, eigenvalues, eigenvectors)

layer=(slice(None),slice(None))

energy=0
orbital=0
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)
fig,ax = plot_data(fig, ax, data).differential_current_map(energy, layer, orbital)
plt.show()
exit()
##########################################################
######################### Plot ###########################
##########################################################
y0=np.real(data.hartree_iterations[0])
y1=np.real(data.gorkov_iterations[0])
print(y0[-1])
print(y1[-1])

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(y0)
ax.plot(y1)
plt.show()
