from lib import *

def main():
    A=Atom([0,0],'A')
    B=Atom([0.5,0],'B')
    A.add_orbital('s')
    B.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    bdg=BogoliubovdeGennes(lattice_vectors,'2D-Weyl-SSH')
    bdg.add_atom(A)
    bdg.add_atom(B)
    bdg.n_spins=2

    bdg.cut(n_cells, axes=0, glue_edgs=True)
    bdg.cut(n_cells, axes=1, glue_edgs=True)
    bdg.set_onsite(-mu+s,atom='A')
    bdg.set_onsite(-mu-s,atom='B')

    bdg.set_hopping(-t,hop_vector=[0,0],atom_i='A',atom_f='B',label='$t$')
    bdg.set_hopping(-t,hop_vector=[-1,0],atom_i='A',atom_f='B',label='$t$')
    bdg.set_hopping(-t,hop_vector=[0,-1],atom_i='A',atom_f='A',label='$t$')
    bdg.set_hopping(-t,hop_vector=[0,-1],atom_i='B',atom_f='B',label='$t$')

    bdg.set_hartree(rho+rho_shift,spin='up')
    bdg.set_hartree(rho-rho_shift,spin='dn')
    bdg.set_fock(phi,atom_i='A',atom_f='B',spin_i='up',spin_f='up')
    bdg.set_fock(phi,atom_i='A',atom_f='B',spin_i='dn',spin_f='dn')
    bdg.set_gorkov(chi+chi_shift,atom_i='A',atom_f='B',spin_i='up',spin_f='up')
    bdg.set_gorkov(chi-chi_shift,atom_i='A',atom_f='B',spin_i='dn',spin_f='dn')

    bdg.set_hubbard_u(Ut,atom_i='A',atom_f='B',hop_vector=[0,0])
    bdg.set_hubbard_u(-Us,atom_i='A',atom_f='A',hop_vector=[0,0],spin_i='up',spin_f='dn')
    bdg.set_hubbard_u(-Us,atom_i='B',atom_f='B',hop_vector=[0,0],spin_i='up',spin_f='dn')

    bdg.record_hartree(location=[0,0], atom='A', spin='up', _print=False)
    bdg.record_hartree(location=[0,0], atom='A', spin='dn', _print=False)
    bdg.record_hartree(location=[0,0], atom='B', spin='up', _print=False)
    bdg.record_hartree(location=[0,0], atom='B', spin='dn', _print=False)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B',spin_i='up', spin_f='up', _print=False)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B',spin_i='dn', spin_f='dn', _print=False)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='A',spin_i='up', spin_f='dn', _print=False)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], atom_i='B', atom_f='B',spin_i='up', spin_f='dn', _print=False)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B',spin_i='up',spin_f='up',_print=False)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='B',spin_i='dn',spin_f='dn',_print=False)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='A', atom_f='A',spin_i='up',spin_f='dn',_print=False)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='B', atom_f='B',spin_i='up',spin_f='dn',_print=False)
        
    if print_free_energy:
        bdg.print_V=True
        bdg.print_V_mf=True
        bdg.print_Eg=True
        bdg.print_free_energy=True

    return bdg

def plot_iterations(bdg):

    markers=['o','+','^','x','.']
    s=3

    fig, ax1  = plt.subplots(1,1,sharex='col')

    # color = 'tab:black'
    ax1.tick_params(axis='y')
    # ax1.plot(bdg._fock_iterations[0],c='k',marker=markers[2],markersize=s,label=f'$\Phi$')
    ax1.plot(bdg._gorkov_iterations[0],c='cyan',marker=markers[3],markersize=s,label=r'$\Delta^\text{AB}_{\uparrow\uparrow}$')
    ax1.plot(bdg._gorkov_iterations[1],c='k',marker=markers[3],markersize=s,label=r'$\Delta^\text{AB}_{\downarrow\downarrow}$',linestyle='dashed')
    ax1.plot(bdg._hartree_iterations[0],c='b',marker=markers[0],markersize=s,label=r'$\phi^\text{A}_\uparrow$')
    ax1.plot(bdg._hartree_iterations[1],c='g',marker=markers[1],markersize=s,label=r'$\phi^\text{A}_\downarrow$')
    ax1.plot(bdg._hartree_iterations[0],c='y',marker=markers[0],markersize=s,label=r'$\phi^\text{B}_\uparrow$',linestyle='dashed')
    ax1.plot(bdg._hartree_iterations[1],c='orange',marker=markers[1],markersize=s,label=r'$\phi^\text{B}_\downarrow$',linestyle='dashed')
    ax1.legend()

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('Free energy', color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.plot(bdg.free_energy,c='r',marker=markers[4],markersize=s,label=f'Free energy')
    # ax2.plot(bdg.Eg,c='b',marker=markers[4],markersize=s,label=f'Eg')
    # ax2.plot(bdg.V,c='g',marker=markers[4],markersize=s,label=f'V')
    # ax2.plot(bdg.V_mf,c='c',marker=markers[4],markersize=s,label=r'V_{mf}')

    # plt.plot(np.real(bdg._gorkov_iterations[1]))
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    xlabel=r'Iterations'
    ylabel=r'Amplitude of fields'
    title=r'BCS theory'
    fig.suptitle(title)
    fig.supxlabel(xlabel)
    fig.supylabel(ylabel)
    plt.tight_layout()
    
#############################################################################
################################# Main ######################################
#############################################################################
mu=4.8
s=0
t=1
Ut=0#2.1
Us=3.2
rho=-2.51
rho_shift=0.179
phi=0#.3
chi=0#3.4
chi_shift=0
n_cells=11

friction=0.
absolute_convergence_factor=0.00001
max_iterations=100
print_free_energy=False

_print=False
bdg = main()
bdg.self_consistent_calculation(friction=friction, max_iterations=400, absolute_convergence_factor=absolute_convergence_factor)

hartree_A_up=bdg.hartree(atom='A',spin='up')
hartree_A_dn=bdg.hartree(atom='A',spin='dn')
hartree_B_up=bdg.hartree(atom='B',spin='up')
hartree_B_dn=bdg.hartree(atom='B',spin='dn')
M_A=(hartree_A_up-hartree_A_dn)[0,0]
M_B=(hartree_B_up-hartree_B_dn)[0,0]
print(M_A)
print(M_B)
plot_iterations(bdg)
plt.show()
exit()

_print=True
print_free_energy=True

# for alpha in [True, False]:
    # greens_function_xy, greens_function_xq, greens_function_kq, bdg = main()
    # greens_function_kq = main()

    # k-space
    # fig,ax = plt.subplots(1,1)
    # greens_function_kq.plot_ldos(ax, energy=0, anomalous=False)

    # iterations
    # fig,ax = plt.subplots(1,1)
# plt.show()

bulk_calculation = False
bdg = main()
# exit()
plot_iterations(bdg)
plt.show()
