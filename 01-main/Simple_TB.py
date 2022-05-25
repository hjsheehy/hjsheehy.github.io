from lib import *

def main():
    A=Atom([0,0],'A')
    A.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    tb=TightBinding(lattice_vectors,'Simple TB')
    tb.add_atom(A)
    tb.n_spins=1
    
    if alpha:
        tb.set_kpts([n_cells,n_cells])
    else:
        tb.cut(n_cells, axes=0, glue_edgs=True)
        tb.cut(n_cells, axes=1, glue_edgs=True)
    tb.set_onsite(-mu)
    tb.set_hopping(-t*1.1,hop_vector=[1,0],label='$t$')
    tb.set_hopping(-t*0.9,hop_vector=[0,1],label='$t$')

    # tb.add_impurities(V,[0,0])

    tb.solve()
    
    energy_interval=np.linspace(-4,4,601)
    resolution=0.05

    # greens_function_xy=GreensFunction(tb,energy_interval,resolution, k_axes=None)

    # greens_function_xq=GreensFunction(tb,energy_interval,resolution, k_axes=[1])

    greens_function_kq=GreensFunction(tb,energy_interval,resolution, k_axes=[0,1])

    del tb.eigenvectors
    del tb.eigenvalues

    return greens_function_kq
    # return greens_function_xy, greens_function_xq, greens_function_kq

# Plotting:

def plot_iterations(bdg):

    markers=['o','+','^','x','.']
    s=3

    # hartree_A=bdg.hartree(atom='A')[0,0]
    # hartree_B=bdg.hartree(atom='B')[0,0]
    # fock_v=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])[0,0]
    # gorkov_v=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[0,0])[0,0]
    # free_energy=bdg.free_energy
    
    # exit()
    # gorkov_w=bdg.gorkov(atom_i='A', atom_f='B', hop_vector=[-1,0])

    fig, [ax1, ax3] = plt.subplots(2,1,sharex='col')

    # color = 'tab:black'
    ax1.tick_params(axis='y')
    ax1.plot(bdg._fock_iterations[0],c='k',marker=markers[2],markersize=s,label=f'$\Phi$')
    ax1.plot(bdg._gorkov_iterations[0],c='cyan',marker=markers[3],markersize=s,label=f'$\Delta$')
    ax1.legend()

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('Free energy', color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.plot(bdg.free_energy,c='r',marker=markers[4],markersize=s,label=f'Free energy')
    # ax2.plot(bdg.Eg,c='b',marker=markers[4],markersize=s,label=f'Eg')
    # ax2.plot(bdg.V,c='g',marker=markers[4],markersize=s,label=f'V')
    # ax2.plot(bdg.V_mf,c='c',marker=markers[4],markersize=s,label=r'V_{mf}')

    ax3.plot(bdg._hartree_iterations[0],c='b',marker=markers[0],markersize=s,label=r'$\phi_\uparrow$')
    ax3.plot(bdg._hartree_iterations[1],c='g',marker=markers[1],markersize=s,label=r'$\phi_\downarrow$')
    ax3.legend()
    # plt.plot(np.real(bdg._gorkov_iterations[1]))
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    xlabel=r'Iterations'
    ylabel=r'Amplitude of fields'
    title=r'Stoner theory'
    fig.suptitle(title)
    fig.supxlabel(xlabel)
    fig.supylabel(ylabel)
    plt.tight_layout()
    
#############################################################################
################################# Main ######################################
#############################################################################
mu=0
t=1
V=0
n_cells=21

bdg = main()

# i=0
# for alpha in [True,False]:
    # greens_function_xy, greens_function_xq, greens_function_kq = main()
    # greens_function_kq = main()

    #real space
    # fig,ax = plt.subplots(1,1)
    # greens_function_xy.plot_ldos(ax, energy=0)
    # plt.show()

    # k-space
    # greens_function_kq.plot_ldos(ax[i], energy=0)
    # i+=1

plot_iterations(bdg)
plt.show()
