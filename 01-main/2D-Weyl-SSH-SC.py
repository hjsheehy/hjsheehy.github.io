from lib import *

filename=sys.argv[0].split('.')[0]
DATA=os.path.join(DATA,filename+'.npz')
FIG=os.path.join(FIG,filename)

def model(*args):
    A=Atom([0,0],'A')
    B=Atom([0.5,0],'B')
    A.add_orbital('s')
    B.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    bdg=BogoliubovdeGennes(lattice_vectors,'2D-Weyl-SSH')
    bdg.add_atom(A)
    bdg.add_atom(B)
    bdg.n_spins=1

    bdg.cut(n_cells, axes=0, glue_edgs=False)
    bdg.cut(n_cells, axes=1, glue_edgs=True)
    bdg.set_onsite(-mu+s,atom='A')
    bdg.set_onsite(-mu-s,atom='B')

    bdg.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B',label='$v$')
    bdg.set_hopping(-w,hop_vector=[-1,0],atom_i='A',atom_f='B',label='$w$')
    bdg.set_hopping(-td,hop_vector=[0,1],atom_i='B',atom_f='A',label='$t_d$')
    bdg.set_hopping(-td,hop_vector=[0,-1],atom_i='B',atom_f='A',label='$t_d$')
    bdg.set_hopping(-td,hop_vector=[1,1],atom_i='B',atom_f='A',label='$t_d$')
    bdg.set_hopping(-td,hop_vector=[1,-1],atom_i='B',atom_f='A',label='$t_d$')

    # impurity_wall = [[0,i] for i in range(n_cells)]
    # bdg.add_impurities(V,impurity_wall)

    # bdg.add_impurities(V,[0,0])

    bdg.set_hartree(rho+rho_shift,atom='A')
    bdg.set_hartree(rho-rho_shift,atom='B')
    bdg.set_fock(phi_v,atom_i='A',atom_f='B')
    bdg.set_fock(phi_w,atom_i='A',atom_f='B',hop_vector=[-1,0])
    bdg.set_gorkov(chi_v,atom_i='A',atom_f='B')
    bdg.set_gorkov(chi_w,atom_i='A',atom_f='B',hop_vector=[-1,0])

    bdg.set_hubbard_u(-Uv,atom_i='A',atom_f='B',hop_vector=[0,0])
    bdg.set_hubbard_u(Uw,atom_i='A',atom_f='B',hop_vector=[-1,0])

    bdg.record_hartree(location=[0,0], atom='A', _print=False)
    bdg.record_hartree(location=[0,0], atom='B', _print=False)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B',_print=False)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B',_print=False)
    # bdg.record_gorkov(location_i=[0,0], location_f=[1,0], atom_i='B', atom_f='A', orbital_i='s', orbital_f='s', spin_i=0, spin_f=0,_print=False)

    return bdg

def process():

    bdg = model()

    bdg.self_consistent_calculation(friction=0.2, max_iterations=200, absolute_convergence_factor=0.00001)

    energy_interval=np.linspace(-4,4,601)
    resolution=0.05

    greens_function_xy=GreensFunction(bdg,energy_interval,resolution, k_axes=None)

    greens_function_xq=GreensFunction(bdg,energy_interval,resolution, k_axes=[1])

    greens_function_kq=GreensFunction(bdg,energy_interval,resolution, k_axes=[0,1])

    del bdg.eigenvectors
    del bdg.eigenvalues

    with open(DATA, 'wb') as f:
        cPickle.dump([greens_function_xy, greens_function_xq, greens_function_kq, bdg], f)
    return greens_function_xy, greens_function_xq, greens_function_kq, bdg

def self_consistent(init_bdg,*args):
    
    bdg = model(*args)
    bdg._hartree=init_bdg._hartree
    bdg._fock=init_bdg._fock
    bdg._gorkov=init_bdg._gorkov
    bdg._hubbard_indices=init_bdg._hubbard_indices
    bdg._anomalous_indices=init_bdg._anomalous_indices
    bdg.U_entries=init_bdg.U_entries
    
    bdg.self_consistent_calculation(friction=0.9, max_iterations=400, absolute_convergence_factor=0.00001)

    hartree_A=bdg.hartree(atom='A')[0,0]
    hartree_B=bdg.hartree(atom='B')[0,0]
    fock_v=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])[0,0]
    gorkov_v=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])[0,0]
    
    del bdg.eigenvectors
    del bdg.eigenvalues

    return bdg, [hartree_A, hartree_B, fock_v, gorkov_v]

def phase_diagrams(xx,zz):
    for z in zz:
        phase_diagram(xx,z)

def phase_diagram(xx,z,include_reverse=True):

    self.init_friction=self.friction
    self.init_max_iterations=self.max_iterations
    self.init_absolute_convergence_factor=self.absolute_convergence_factor

    self.init_friction=0.7
    self.init_max_iterations=400
    self.init_absolute_convergence=0.00001
    
    yy=[]
    yyy=[]
    if include_reverse:
        xxx=[xx, xx[::-1]]
    else:
        xxx=[xx]

    for xx in xxx:
        for i,x in enumerate(xx):
            if i==0:
                bdg = model(x,z)
                bdg.self_consistent_calculation(friction=self.init_friction, max_iterations=self.init_max_iterations, absolute_convergence_factor=self.init_absolute_convergence_factor)
            bdg, y = self_consistent(bdg,x,z)
            yy.append(y)
        yyy.append(yy)

    xxx=np.array(xxx)
    yyy=np.array(yyy)

    with open(DATA+'_phase_diagram_{Uv:.2f}', 'wb') as f:
        cPickle.dump([xxx,yyy,z], f)

def plot_phase_diagram():
    [Uv,fields] = np.load(DATA+'_phase_diagram', allow_pickle=True)
    

    print(np.shape(fields))
    markers=['>','<']
    plt.plot(muu,fields[0,:,3],marker=markers[0])
    plt.plot(muu,fields[1,:,3],marker=markers[1])
    plt.show()


# Plotting:

def plot_iterations(bdg):
    hartree_A=bdg.hartree(atom='A')[0,0]
    hartree_B=bdg.hartree(atom='B')[0,0]
    fock_v=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])[0,0]
    gorkov_v=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])[0,0]
    
    # exit()
    # gorkov_w=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[-1,0])

    plt.plot(bdg._hartree_iterations[0],c='b')
    plt.plot(bdg._hartree_iterations[1],c='g')
    plt.plot(bdg._fock_iterations[0],c='k')
    plt.plot(bdg._gorkov_iterations[0],c='r')
    # plt.plot(np.real(bdg._gorkov_iterations[1]))
    plt.show()
    plt.close()

    # FIGNAME='unit_cell'

    # fig, ax = plt.subplots(1, 1)
    # fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    # plt.tight_layout()
    # model.plot_unit_cell(fig, ax, atoms='all', s=100)
    # output=os.path.join(FIG,FIGNAME)
    # plt.savefig(output+'.pdf', bbox_inches = "tight")

    # with open(output+'.txt', 'w') as f:
    #     f.write(rf'''The local density of states of a spinless square lattice tight-binding
# model at $\mu/t={mu:.2f}$ with ${nx}\times{ny}$
# sites, a single orbital with an impurity at the centre with coupling
# strength $V/t={V:.2f}$. The impurity gives rise to Fridel's eponymous
# waves in the electron quasiparticle density. The electronic excitations
# at zero temperature necessarily carry the Fermi energy, and hence the 
# wavefunction describing the excitations is of the Fermi wavelength. 
# The electronic charge distrbution is the square modulus of the
# wavefunction and hence takes on twice the periodicty or double the
# wavelength $\lambda_\text{{Friedel}}=\lambda_\text{{Fermi}}/2=
# {friedel_wavelength:.3}$. \\
# ''')

#############################################################################
################################# Main ######################################
#############################################################################
mu=2.7
s=0.0
td=0.9
v=0.6
w=1.2
Uv=3.6
Uw=0
rho=6.0
rho_shift=0.3
phi_v=0.2
phi_w=0
chi_v=4.2
chi_w=0
V=0
n_cells=3

# greens_function_xy, greens_function_xq, greens_function_kq, bdg = process()

# [greens_function_xy, greens_function_xq, greens_function_kq, bdg] = np.load(DATA, allow_pickle=True)

# plot_iterations(bdg)

mu=np.arange(-4,4.1,0.1)
Uv=np.arange(1.1,6.5,1.1)
phase_diagrams(mu,Uv)

plot_phase_diagram()
