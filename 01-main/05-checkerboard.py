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
temperature=0
max_iterations=500
friction=0.8
absolute_convergence_factor=0.0001
############## Energy ###############
resolution=0.1
increment=0.01
max_val=10
energy_interval = np.arange(start=-max_val, stop=max_val+increment, step=increment)
##########################################################
######################### Main ###########################
##########################################################
model = bogoliubov_de_gennes(dimensions=[n,1], n_spins=2, basis=[[1,0],[0,1]], orbitals=[[0.35,0.5],[0.65,0.5]], pbc=[True,True])
model.set_onsite(-mu)
model.set_hopping(-v, 0, 1, [0,0],'$-v$')
model.set_hopping(-w, 1, 0, [1,0],'$-w$')
model.set_hopping(-w, 0, 0, [0,1],'$-w$')
model.set_hopping(-w, 1, 1, [0,1],'$-w$')
model.set_impurities(V, locations=[0],spins=[0,1],orbitals=[0,1])

model.set_hubbard_u(-U_T*np.eye(2), orb_i=0, orb_f=1, hop_vector=[0,0])
hubbard_R=U_R*np.eye(2)
model.set_hubbard_u(hubbard_R, orb_i=1, orb_f=0, hop_vector=[1,0])
model.set_hubbard_u(hubbard_R, orb_i=0, orb_f=0, hop_vector=[0,1])
model.set_hubbard_u(hubbard_R, orb_i=1, orb_f=1, hop_vector=[0,1])

model.set_hartree([1.2*phi,0.8*phi])

model.set_fock(zeta_v*np.eye(2), orb_i=0, orb_f=1, hop_vector=[0,0])
model.set_fock(zeta_w*np.eye(2), orb_i=1, orb_f=0, hop_vector=[1,0])
model.set_fock(zeta_w*np.eye(2), orb_i=0, orb_f=0, hop_vector=[0,1])
model.set_fock(zeta_w*np.eye(2), orb_i=1, orb_f=1, hop_vector=[0,1])

model.set_gorkov(Delta_v*np.eye(2), orb_i=0, orb_f=1, hop_vector=[0,0])
model.set_gorkov(Delta_w*np.eye(2), orb_i=1, orb_f=0, hop_vector=[1,0])
model.set_gorkov(Delta_w*np.eye(2), orb_i=0, orb_f=0, hop_vector=[0,1])
model.set_gorkov(Delta_w*np.eye(2), orb_i=1, orb_f=1, hop_vector=[0,1])

######## Importing previous meanfields #######
dataDir=os.listdir(os.path.join(DATA,fileName))

if len(dataDir)>0:
    MU=np.array([float(txt.split("_")[1]) for txt in dataDir])
    R=np.array([float(txt.split("_")[2]) for txt in dataDir])
    T=np.array([float(txt.split("_")[3]) for txt in dataDir])

    point=[mu,U_R,U_T]
    points=np.array([MU,R,T]).T
    i=closest_point(points,point)
    [MU,R,T]=[MU[i],R[i],T[i]]
    val=f'{mu:0.2f}_{R:0.4f}_{T:0.4f}'
    i=len(dataDir)
    tmp=False
    for ss in dataDir:
        i-=1
        if any(val in s for s in dataDir):
            closestFile=dataDir[i]
            tmp=True
    if tmp:
        closestFile = os.path.join(DATA,fileName,closestFile)
        with open(closestFile, 'rb') as f:
            data = cPickle.load(f)
        
        hartree=data.hartree()
        fock=data.fock()
        gorkov=data.gorkov()

        if np.any(np.isnan(hartree)):
            print('is nan')

            hartree=np.random.uniform(low=0.1, high=1.3, size=np.shape(data.hartree()))
            fock=np.random.uniform(low=0.1, high=2.3, size=np.shape(data.fock()))
            gorkov=np.random.uniform(low=0.1, high=4.3, size=np.shape(data.gorkov()))

            extended_dimensions=np.append(model.extended_dimensions,model.extended_dimensions)

            model.hartree=np.reshape(hartree,model.n_dof,'F')
            model.fock=np.reshape(fock,[model.n_dof,model.n_dof],'F')
            model.gorkov=np.reshape(gorkov,[model.n_dof,model.n_dof],'F')

        extended_dimensions=np.append(model.extended_dimensions,model.extended_dimensions)

        model.hartree=np.real_if_close(np.reshape(hartree,model.n_dof,'F'))
        model.fock=np.real_if_close(np.reshape(fock,[model.n_dof,model.n_dof],'F'))
        model.gorkov=np.real_if_close(np.reshape(gorkov,[model.n_dof,model.n_dof],'F'))

_print=False
# _print=True

# model.record_hartree([0,0], 0, 0, _print)
# model.record_hartree([0,0], 1, 1, _print)
model.record_gorkov(location_a=[5,0], location_b=[5,0], spin_a=0, spin_b=0, orbital_a=0, orbital_b=1, _print=_print)
model.record_gorkov(location_a=[4,0], location_b=[5,0], spin_a=0, spin_b=0, orbital_a=1, orbital_b=0, _print=_print)
# model.record_gorkov(location_a=[0,0], location_b=[0,1], spin_a=0, spin_b=0, orbital_a=0, orbital_b=0, _print=_print)
# model.record_fock(location_a=[5,0], location_b=[5,0], spin_a=0, spin_b=0, orbital_a=0, orbital_b=1, _print=_print)
# model.record_fock(location_a=[4,0], location_b=[5,0], spin_a=0, spin_b=0, orbital_a=1, orbital_b=0, _print=_print)
# model.record_fock(location_a=[0,0], location_b=[1,0], spin_a=0, spin_b=0, orbital_a=1, orbital_b=0, _print=_print)
# model.record_fock(location_a=[0,0], location_b=[1,0], spin_a=0, spin_b=1, orbital_a=1, orbital_b=0, _print=_print)

model.set_max_iterations(max_iterations)
model.set_friction(friction)
model.set_absolute_convergence_factor(absolute_convergence_factor)
model.set_temperature(temperature)
eigenvalues,eigenvectors=model.self_consistent_calculation()    

reset=False
attempts=1
for counter in range(attempts):
    if not model.converged and reset:
        print("Repeating renormalisation with reset fields.")

        model.reset_hartree()
        model.reset_fock()
        model.reset_gorkov()

        hartree=np.random.uniform(low=0.1, high=1.3, size=np.shape(data.hartree()))
        fock=np.random.uniform(low=0.1, high=2.3, size=np.shape(data.fock()))
        gorkov=np.random.uniform(low=0.1, high=4.3, size=np.shape(data.gorkov()))

        extended_dimensions=np.append(model.extended_dimensions,model.extended_dimensions)

        model.hartree=np.reshape(hartree,model.n_dof,'F')
        model.fock=np.reshape(fock,[model.n_dof,model.n_dof],'F')
        model.gorkov=np.reshape(gorkov,[model.n_dof,model.n_dof],'F')

        eigenvalues,eigenvectors=model.self_consistent_calculation()    

nan=True
if not model.converged and nan:
    model.hartree=model.hartree*float('nan')
    model.fock=model.fock*float('nan')
    model.gorkov=model.gorkov*float('nan')
        
# Trash unnecessary data
del(model._tb_ham)
del(model._hubbard_u)

data=processed_data(model, energy_interval, resolution, eigenvalues, eigenvectors)

with open(confname+'.npz', 'wb') as f:
    cPickle.dump(data, f)
exit()
########################################################
########################## Plot ########################
########################################################
latex_width=4.7747
n_pts=40    
def fig_a():
    energy=0

    fig, axs = plt.subplots(1)
    plt_data = plot_data(fig,axs,data)

    fig, axs = plt_data.differential_current_map(energy, layer=(ALL,ALL), orbital=[0,1], spin=None, cartesian=True, n_pts=n_pts, gaussian_mean=0.1)

    label_x = r'$x/a_x$'
    label_y = r'$y/a_y$'
    axs.set_ylabel(label_y)
    axs.set_xlabel(label_x)

    plt.tight_layout()

    axs.set_title('zero-bias local density of states')
    return fig
y=data.centre[1]
gorkov=data.gorkov()
Delta_T=np.real([[gorkov[x,y,0,s,x,y,1,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
Delta_T=np.sum(Delta_T,0)/data.n_spins
Delta_R=np.real([[gorkov[x-1,y,1,s,x,y,0,s] for x in range(data.dimensions[0])] for s in range(data.n_spins)])
Delta_R=np.sum(Delta_R,0)/data.n_spins
print(f'Delta_T={Delta_T[-1]}')
print(f'Delta_R={Delta_R[-1]}')

fig_a()
plt.show()
