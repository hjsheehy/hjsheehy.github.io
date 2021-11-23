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
basis=[[1,0],[0,1]]
pbc=[bc,False]
orbitals=[[0,0],[0.25,0.5]]
n_orbs=len(orbitals)
dimensions=[n_x,1]
orbitals=[[0,0],[0.25,0.5]]
hoppings=[[-t, 0, 0, [1,0],'$-t$'],[-t, 1, 1, [1,0],'$-t$'],[-v, 0, 1, [0,0],'$-v$'],[-w, 1, 0, [1,0], '$-w$']]
mu_a=mu-s
mu_b=mu+s
SO_spin=np.eye(n_spins)
SO_orbital=np.eye(n_orbs)
SO_onsite=-mu*kron(SO_spin,SO_orbital)
impurity_location = [[0,0]]
M=[0,0,1]
############## Energy ###############
resolution=0.1
increment=0.01
max_val=10
energy_interval = np.arange(start=-max_val, stop=max_val+increment, step=increment)
##########################################################
######################### Main ###########################
##########################################################
if model=='tb':
    model = tight_binding(dimensions, n_spins, basis, orbitals, pbc)
    #model.set_SO_onsite(SO_onsite)
    model.set_onsite(-mu_a,orbital=0)
    model.set_onsite(-mu_b,orbital=1)
    for link in hoppings:
        model.set_hopping(*link)
    #model.set_impurities(V, impurity_location)
    model.set_magnetic_impurities(M,impurity_location)
    axes=[0,1]
    #model.fourier_transform_hamiltonian(transform=axes)
    #model.fourier_transform_hamiltonian(inverse_transform=axes)
    eigenvalues,eigenvectors=model.solve()
    data=processed_data(model, energy_interval, resolution, eigenvalues, eigenvectors)

elif model=='mf':
    model = bogoliubov_de_gennes(dimensions, n_spins, basis, orbitals, pbc)
    model.set_onsite(-mu_a,spin=0)
    model.set_onsite(-mu_b,spin=1)
    for link in hoppings:
        model.set_hopping(*link)
    #model.set_impurities(V, impurity_location)
    #model.set_hartree([1.2*phi,0.8*phi])
    spin_tensor=Delta*(1.0j*Pauli_y)
    model.set_gorkov(spin_tensor)
    model.set_mean_field_hamiltonian()
    eigenvalues,eigenvectors=model.solve()
    data=processed_data(model, energy_interval, resolution, eigenvalues, eigenvectors)

with open(confname+'.npz', 'wb') as f:
    cPickle.dump(data, f)
exit()
########################################################
######################## temp ##########################
########################################################
import matplotlib.pyplot as plt
fig, axs = plt.subplots(1)
latex_width=4.7747
fig.set_size_inches(w=latex_width, h=4.5) 

energy=0
layer=(slice(None),slice(None))
plot_data = plot_data(fig,axs,data)
#dos = plot_data.density_of_states
orbital=0

fig,axs = plot_data.lattice(energy=energy, model=data, layer=layer, spins=np.arange(n_spins), orbitals=np.arange(data.n_orbitals), annotate_orbs=True, show_cell_borders=True, show_basis_vectors=False, show_hopping=True, show_impurities=True, show_only_centre=True)

########################################################
########################## Plot ########################
########################################################
def main():
    label_x = r'$x/a_x$'
    label_y = r'$y/a_y$'
    
    i=0
    axs.set_ylabel(label_y)

    axs.set_xlabel(label_x)

    plt.tight_layout()
    return fig
#fig = main()
#energy
dimension=0
#fig,axs = plot_data.band_structure(dimension)

#fig, axs = plot_data.differential_current_map(energy, layer=(ALL,ALL), orbital=[0,1], spin=None, cartesian=True, n_pts=40, gaussian_mean=0.2)
plt.show()
# plt.savefig('test.pdf', bbox_inches = "tight")
