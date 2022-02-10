from lib import *
#########################################################
#################### BCS wire model #####################
#########################################################
dimensions=[1]
mu=-.7
t=1
n_spins=2
orbitals=[[0]]
basis=[[1]]
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
model.set_hopping(-t, 0, 0, [1])
model.set_impurities(0, [0])
model.set_hubbard_u(hubbard_U, None, None, [0])
model.set_hartree([phiUp,phiDown])
model.set_gorkov(Delta*(1.0j*Pauli_y))
_print=False

model.record_hartree([0], 0, 0, _print)
model.record_gorkov([0], [0], 0, 1, 0, 0, _print)

model.set_max_iterations(max_iterations=10000)
model.set_friction(0)
model.set_temperature(0.05)
model.set_absolute_convergence_factor(0.00001)
eigenvalues,eigenvectors=model.self_consistent_calculation()    
data=processed_data(model, energy_interval, resolution, eigenvalues, eigenvectors)
print(model.free_energy[-1])
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
