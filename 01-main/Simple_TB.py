from lib import *

def main():
    A=Atom([0,0],'A')
    A.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    tb=TightBinding(lattice_vectors,'Simple TB')
    tb.add_atom(A)
    tb.n_spins=2
    
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

#############################################################################
################################# Main ######################################
#############################################################################
mu=-3.2
t=1
V=0
n_cells=21

fig,ax = plt.subplots(1,2)
i=0
for alpha in [True,False]:
    # greens_function_xy, greens_function_xq, greens_function_kq = main()
    greens_function_kq = main()

    #real space
    # fig,ax = plt.subplots(1,1)
    # greens_function_xy.plot_ldos(ax, energy=0)
    # plt.show()

    # k-space
    greens_function_kq.plot_ldos(ax[i], energy=0)
    i+=1
plt.show()
