from PhaseDiagram import *

def model(U,L):

    Cu=Atom([0,0],'Cu^{d}')
    O_x=Atom([0.5,0],'O^{p_x}')
    O_y=Atom([0,0.5],'O^{p_y}')
    Cu.add_orbital('d')
    O_x.add_orbital('p_x')
    O_y.add_orbital('p_y')
    lattice_vectors=[[1,0],[0,1]]
    bdg=BogoliubovdeGennes(lattice_vectors,'CuO')
    bdg.add_atom(A)
    bdg.n_spins=2

    bdg.cut(n_cells, [0,1], glue_edgs=True)
    bdg.set_onsite(-mu)
    bdg.set_hopping(-t,atom_i='Cu',atom_f='O_x',hop_vector=[0,0],label='t')
    bdg.set_hopping(-t,atom_i='Cu',atom_f='O_y',hop_vector=[0,0],label='t')
    bdg.set_hopping(-t,atom_i='O_x',atom_f='Cu',hop_vector=[1,0],label='t')
    bdg.set_hopping(-t,atom_i='O_y',atom_f='Cu',hop_vector=[0,1],label='t')

    impurity_locations=points_in_circle(radius=L, x0=0, y0=0)
    
    bdg.add_impurities(V,impurity_locations,label='V')
    
    if antiferromagnetic:
        bdg.set_hartree_antiferromagnetic(rho,rho_shift)
    else:
        bdg.set_hartree(rho+rho_shift,spin='up')
        bdg.set_hartree(rho-rho_shift,spin='dn')
    bdg.set_fock(phi,spin_i='up',spin_f='dn')
    bdg.set_gorkov(chi,spin_i='up',spin_f='dn')
    
    _print=False
    ##############################
    ### Add to notes and docs! ###
    ###     U>0 repulsive      ###
    ###     U<0 attractive     ###
    ##############################
    bdg.set_hubbard_u(U,spin_i='up',spin_f='dn',hop_vector=[0,0])
    bdg.set_hubbard_u_impurities(-U, atom='Cu', impurity_locations=impurity_locations)
    bdg.U=U

    bdg.record_hartree(location=[0,0], spin='up', _print=_print)
    bdg.record_hartree(location=[0,0], spin='dn', _print=_print)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='dn',_print=_print)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='dn',_print=_print)

    return bdg

def func_of_dep_vars(bdg):
    energy_interval=np.linspace(-8,8,600)
    resolution=0.1
    bdg.greens_function_kq=GreensFunction(bdg,energy_interval,resolution, k_axes=[0,1])

    return bdg

def plotting_func(bdg):
    fig,ax=plt.subplots()
    # ax = bdg.greens_function_kq.plot_spectrum(ax, energy='resolved', axes=['resolved','integrated'], atom='integrated', xmin='default', xmax='default', omega_min='default', omega_max='default', vmin=0, vmax=60,label='')
    
    # Plot iterations:
    # plt.plot(np.real(bdg._hartree_iterations[0]),color='b',label=r'$\phi_\uparrow$')
    # plt.plot(np.real(bdg._hartree_iterations[1]),color='r',label=r'$\phi_\downarrow$')
    # plt.plot(np.real(bdg._gorkov_iterations[0]),color='g',label=r'$\Delta_{\uparrow\downarrow}$')
    # plt.legend()

    # plot DOS:
    greens_function = bdg.greens_function_kq

    x = greens_function.energy_interval
    y = greens_function.spectrum(sites='integrated', atom='integrated',orbital='integrated', spin='integrated', anomalous=False)
    ax.plot(y,x)
    ax.set_xlabel(r'$\text{DOS}(\omega)$')
    ax.set_ylabel(r'$\omega$')
    N=(bdg.hartree(spin='up')+bdg.hartree(spin='dn'))[0,0]
    M=np.abs(bdg.hartree(spin='up')-bdg.hartree(spin='dn'))[0,0]
    Delta=bdg.gorkov(spin_i='up',spin_f='dn')[0,0]
    ax.set_title(rf'''$U={bdg.U:0.03f}$
$N={N:0.03f}, M={M:0.03f}, \Delta={Delta:0.03f}$''')

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    plt.tight_layout()
    return fig

def minimised_plot_func():
    pass

####################################################

n_cells=33
mu=0
t=1

V=0

U=-4.1
antiferromagnetic=False


bdg = model(*args)

bdg.self_consistent_calculation(friction=friction, max_iterations=400, absolute_convergence_factor=absolute_convergence_factor)

energy_interval=np.linspace(-4,4,601)
resolution=0.05

# greens_function_xy=GreensFunction(bdg,energy_interval,resolution, k_axes=None)

# greens_function_xq=GreensFunction(bdg,energy_interval,resolution, k_axes=[1])

greens_function_kq=GreensFunction(bdg,energy_interval,resolution, k_axes=[0,1])






######################################################

# pos=np.linspace(0,10,20)
# neg=np.linspace(-10,0,20)
# XX=[pos,pos[::-1],neg,neg[::-1]]

# rho,rho_shift,phi,chi=5,0,0,4
# names=['Delta_rep_U_fwd','Delta_rep_U_rev','Delta_att_U_fwd','Delta_att_U_rev']
# for i in range(len(names)):
#     name=names[i]
#     X=XX[i]

#     phase_diagram=PhaseDiagram(model=model)
#     phase_diagram.set_directory_name(name)
#     phase_diagram.set_indep_vars(X)
#     phase_diagram.set_dep_vars_func(func_of_dep_vars)
#     phase_diagram.set_plotting_func(plotting_func)
#     # phase_diagram.execute()
#     phase_diagram.set_png()
#     phase_diagram.plot()

# rho,rho_shift,phi,chi=-2,5,0,0
# names=['Mag_rep_U_fwd','Mag_rep_U_rev','Mag_att_U_fwd','Mag_att_U_rev']
# for i in range(len(names)):
#     name=names[i]
#     X=XX[i]

#     phase_diagram=PhaseDiagram(model=model)
#     phase_diagram.set_directory_name(name)
#     phase_diagram.set_indep_vars(X)
#     phase_diagram.set_dep_vars_func(func_of_dep_vars)
#     phase_diagram.set_plotting_func(plotting_func)
#     # phase_diagram.execute()
#     phase_diagram.set_png()
#     phase_diagram.plot()

# names=['AF_rep_U_fwd','AF_rep_U_rev','AF_att_U_fwd','AF_att_U_rev']
# antiferromagnetic=True
# for i in range(len(names)):
#     name=names[i]
#     X=XX[i]

#     phase_diagram=PhaseDiagram(model=model)
#     phase_diagram.set_directory_name(name)
#     phase_diagram.set_indep_vars(X)
#     phase_diagram.set_dep_vars_func(func_of_dep_vars)
#     phase_diagram.set_plotting_func(plotting_func)
#     # phase_diagram.execute()
#     phase_diagram.set_png()
#     phase_diagram.plot()

# names=['Delta_rep_U_fwd','Delta_rep_U_rev','Delta_att_U_fwd','Delta_att_U_rev','Mag_rep_U_fwd','Mag_rep_U_rev','Mag_att_U_fwd','Mag_att_U_rev','AF_rep_U_fwd','AF_rep_U_rev','AF_att_U_fwd','AF_att_U_rev']
# minimised_data = PlotMinimisedData()
# minimised_data.set_data_directories(names)
# minimised_data.set_directory_name('minimised_data')
# minimised_data.set_plotting_func(plotting_func)
# minimised_data.minimised_data()
# minimised_data.set_png()
# minimised_data.plot()
