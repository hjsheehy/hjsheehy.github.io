from lib import *

def model():

    A=Atom([0,0],'up')
    A.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    bdg=BogoliubovdeGennes(lattice_vectors,'SSH')
    bdg.add_atom(A)
    bdg.n_spins=2

    bdg.cut(n_cells, [0,1], glue_edgs=True)
    bdg.set_onsite(-mu,orbital='s')
    bdg.set_hopping(-t,hop_vector=[1,0],label='t')
    bdg.set_hopping(-t,hop_vector=[0,1],label='t')
    
    bdg.add_impurities(V,impurity_locations,label='V')

    bdg.set_hartree_antiferromagnetic(rho,rho_shift)
    bdg.set_fock(phi,spin_i='up',spin_f='dn')
    bdg.set_gorkov(chi,spin_i='up',spin_f='dn')
    
    _print=False
    bdg.set_hubbard_u(U,spin_i='up',spin_f='dn',hop_vector=[0,0])
    bdg.U=U

    bdg.record_hartree(location=[0,0], spin='up', _print=_print)
    bdg.record_hartree(location=[0,0], spin='dn', _print=_print)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='dn',_print=_print)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='dn',_print=_print)

    return bdg

n_cells=11
mu=0
t=1
V=0
impurity_locations=[[0,0]]
U=3
rho,rho_shift,phi,chi=0,1,0,0
bdg=model()

energy_interval=np.linspace(-5,5,200)
resolution=0.2
greens_function_xy=GreensFunction(bdg,energy_interval,resolution, k_axes=None)

mag=greens_function_xy.magnetism(energy=0)

mag=np.fft.fftshift(mag)
plt.imshow(mag)
plt.colorbar()
plt.show()




