from PhaseDiagram import *

def model(U):

    Cu=Atom([0,0],'Cu')
    O_x=Atom([0.5,0],'O_x')
    O_y=Atom([0,0.5],'O_y')
    Cu.add_orbital('d')
    O_x.add_orbital('p_x')
    O_y.add_orbital('p_y')
    lattice_vectors=[[1,0],[0,1]]
    bdg=BogoliubovdeGennes(lattice_vectors,'CuO_2')
    bdg.add_atom(Cu)
    bdg.add_atom(O_x)
    bdg.add_atom(O_y)
    bdg.n_spins=2

    bdg.cut(n_cells, [0,1], glue_edgs=True)
    bdg.set_onsite(-mu+s,atom='Cu')
    bdg.set_onsite(-mu-s,atom='O_x')
    bdg.set_onsite(-mu-s,atom='O_y')
    bdg.set_hopping(-t,atom_i='Cu',atom_f='O_x',hop_vector=[0,0],label='-t')
    bdg.set_hopping(t,atom_i='Cu',atom_f='O_y',hop_vector=[0,0],label='t')
    bdg.set_hopping(t,atom_i='O_x',atom_f='Cu',hop_vector=[1,0],label='t')
    bdg.set_hopping(-t,atom_i='O_y',atom_f='Cu',hop_vector=[0,1],label='-t')

    impurity_locations=points_in_circle(radius=L, x0=0, y0=0)

    bdg.add_impurities(V,impurity_locations,label='V')

    if antiferromagnetic:
        bdg.set_hartree_antiferromagnetic(rho,rho_shift,atom='Cu')
        bdg.set_hartree_antiferromagnetic(rho,-rho_shift,atom='O_x')
        bdg.set_hartree_antiferromagnetic(rho,-rho_shift,atom='O_y')
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

    # bdg.record_hartree(location=[0,0], spin='up', _print=_print)
    # bdg.record_hartree(location=[0,0], spin='dn', _print=_print)
    # bdg.record_fock(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='dn',_print=_print)
    # bdg.record_gorkov(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='dn',_print=_print)

    return bdg

def func_of_dep_vars(bdg):
    energy_interval=np.linspace(-8,8,600)
    resolution=0.1
    bdg.greens_function_kq=GreensFunction(bdg,energy_interval,resolution, k_axes=[0,1])
    bdg.greens_function_xy=GreensFunction(bdg,energy_interval,resolution, k_axes=None)
    bdg.greens_function_xq=GreensFunction(bdg,energy_interval,resolution, k_axes=[1])

    return bdg

def plotting_func(bdg):
    fig,[ax1,ax2]=plt.subplots(1,2)
    # ax = bdg.greens_function_kq.plot_spectrum(ax, energy='resolved', axes=['resolved','integrated'], atom='integrated', xmin='default', xmax='default', omega_min='default', omega_max='default', vmin=0, vmax=60,label='')

    # Plot iterations:
    # plt.plot(np.real(bdg._hartree_iterations[0]),color='b',label=r'$\phi_\uparrow$')
    # plt.plot(np.real(bdg._hartree_iterations[1]),color='r',label=r'$\phi_\downarrow$')
    # plt.plot(np.real(bdg._gorkov_iterations[0]),color='g',label=r'$\Delta_{\uparrow\downarrow}$')
    # plt.legend()

    # plot DOS:
    greens_function = bdg.greens_function_kq

    x = greens_function.energy_interval
    y1 = greens_function.spectrum(sites='integrated', atom='Cu',orbital='integrated', spin='integrated', anomalous=False)
    y2 = greens_function.spectrum(sites='integrated', atom='O_x',orbital='integrated', spin='integrated', anomalous=False)
    y3 = greens_function.spectrum(sites='integrated', atom='O_y',orbital='integrated', spin='integrated', anomalous=False)
    linestyles=['solid','solid','dashed']
    ax1.plot(y1,x,label=r'$\text{Cu}^{d_{x^2-y^2}}$',c='r',linestyle=linestyles[0])
    ax1.plot(y2,x,label=r'$\text{O}^{p_x}$',c='g',linestyle=linestyles[1])
    ax1.plot(y3,x,label=r'$\text{O}^{p_y}$',c='b',linestyle=linestyles[2])
    ax1.set_xlabel(r'$\text{DOS}(\omega)$')
    ax1.set_ylabel(r'$\omega$')
    ax1.legend()
    N1=(bdg.hartree(atom='Cu',spin='up')+bdg.hartree(atom='Cu',spin='dn'))[0]
    N2=(bdg.hartree(atom='O_x',spin='up')+bdg.hartree(atom='O_x',spin='dn'))[0]
    N3=(bdg.hartree(atom='O_y',spin='up')+bdg.hartree(atom='O_y',spin='dn'))[0]
    M1=np.abs(bdg.hartree(atom='Cu',spin='up')-bdg.hartree(atom='Cu',spin='dn'))[0]
    M2=np.abs(bdg.hartree(atom='O_x',spin='up')-bdg.hartree(atom='O_x',spin='dn'))[0]
    M3=np.abs(bdg.hartree(atom='O_y',spin='up')-bdg.hartree(atom='O_y',spin='dn'))[0]
    Phi1=bdg.fock(atom_i='Cu',atom_f='Cu',spin_i='up',spin_f='dn')[0]
    Phi2=bdg.fock(atom_i='O_x',atom_f='O_x',spin_i='up',spin_f='dn')[0]
    Phi3=bdg.fock(atom_i='O_y',atom_f='O_y',spin_i='up',spin_f='dn')[0]
    Delta1=bdg.gorkov(atom_i='Cu',atom_f='Cu',spin_i='up',spin_f='dn')[0]
    Delta2=bdg.gorkov(atom_i='O_x',atom_f='O_x',spin_i='up',spin_f='dn')[0]
    Delta3=bdg.gorkov(atom_i='O_y',atom_f='O_y',spin_i='up',spin_f='dn')[0]
    x=np.arange(-int(n_cells/2),int(n_cells/2)+1,1)
    markers=['^','v','.','x']
    from matplotlib.lines import Line2D
    lines = [Line2D([0], [0], color='k', marker=marker, linewidth=1, linestyle='solid') for marker in markers]
    labels = [r'$\langle\hat{\mathcal{N}}\rangle$',r'$\langle\hat{\mathcal{M}}\rangle$',r'$\Phi$',r'$\Delta$']
    ax2.legend(lines, labels)
    i=0
    ax2.plot(x,N1,label=r'$\text{Cu}^{d_{x^2-y^2}}$',c='r',linestyle=linestyles[0],marker=markers[i])
    ax2.plot(x,N2,label=r'$\text{O}^{p_x}$',c='g',linestyle=linestyles[1],marker=markers[i])
    ax2.plot(x,N3,label=r'$\text{O}^{p_y}$',c='b',linestyle=linestyles[2],marker=markers[i])
    i=1
    ax2.plot(x,M1,label=r'$\text{Cu}^{d_{x^2-y^2}}$',c='r',linestyle=linestyles[0],marker=markers[i])
    ax2.plot(x,M2,label=r'$\text{O}^{p_x}$',c='g',linestyle=linestyles[1],marker=markers[i])
    ax2.plot(x,M3,label=r'$\text{O}^{p_y}$',c='b',linestyle=linestyles[2],marker=markers[i])
    i=2
    ax2.plot(x,Phi1,label=r'$\text{Cu}^{d_{x^2-y^2}}$',c='r',linestyle=linestyles[0],marker=markers[i])
    ax2.plot(x,Phi2,label=r'$\text{O}^{p_x}$',c='g',linestyle=linestyles[1],marker=markers[i])
    ax2.plot(x,Phi3,label=r'$\text{O}^{p_y}$',c='b',linestyle=linestyles[2],marker=markers[i])
    i=3
    ax2.plot(x,Delta1,label=r'$\text{Cu}^{d_{x^2-y^2}}$',c='r',linestyle=linestyles[0],marker=markers[i])
    ax2.plot(x,Delta2,label=r'$\text{O}^{p_x}$',c='g',linestyle=linestyles[1],marker=markers[i])
    ax2.plot(x,Delta3,label=r'$\text{O}^{p_y}$',c='b',linestyle=linestyles[2],marker=markers[i])
    # ax.set_title(rf'''$U={bdg.U:0.03f}$
# $N={N:0.03f}, M={M:0.03f}, \Delta={Delta:0.03f}$''')

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    plt.tight_layout()
    return fig

def minimised_plot_func():
    pass

####################################################

n_cells=5
mu=0
s=0
t=1

V=0

L=0

rho,rho_shift,phi,chi=10,0,0,3
antiferromagnetic=False
U=4

DATA=DATA+'temp'
# bdg=model(U)
# bdg.self_consistent_calculation(friction=0.9, max_iterations=400, absolute_convergence_factor=0.00001)
# func_of_dep_vars(bdg)
# with open(DATA+'.npz', 'wb') as f:
#     cPickle.dump(bdg, f)
bdg = np.load(os.path.join(DATA+'.npz'), allow_pickle=True)

plotting_func(bdg)
plt.show()

exit()

EXECUTE=True

pos=np.linspace(0,10,20)
neg=np.linspace(-10,0,20)
XX=[pos,pos[::-1],neg,neg[::-1]]

rho,rho_shift,phi,chi=10,0,0,3
antiferromagnetic=False
names=['Delta_rep_U_fwd','Delta_rep_U_rev','Delta_att_U_fwd','Delta_att_U_rev']
for i in range(len(names)):
    name=names[i]
    X=XX[i]

    phase_diagram=PhaseDiagram(model=model)
    phase_diagram.set_directory_name(name)
    phase_diagram.set_indep_vars(X)
    phase_diagram.set_dep_vars_func(func_of_dep_vars)
    phase_diagram.set_plotting_func(plotting_func)
    phase_diagram.set_png()
    if EXECUTE:
        phase_diagram.execute()
    else:
        phase_diagram.plot()

rho,rho_shift,phi,chi=2,10,0,0
antiferromagnetic=False
names=['Mag_rep_U_fwd','Mag_rep_U_rev','Mag_att_U_fwd','Mag_att_U_rev']
for i in range(len(names)):
    name=names[i]
    X=XX[i]

    phase_diagram=PhaseDiagram(model=model)
    phase_diagram.set_directory_name(name)
    phase_diagram.set_indep_vars(X)
    phase_diagram.set_dep_vars_func(func_of_dep_vars)
    phase_diagram.set_plotting_func(plotting_func)
    phase_diagram.set_png()
    if EXECUTE:
        phase_diagram.execute()
    else:
        phase_diagram.plot()

antiferromagnetic=True
names=['AF_rep_U_fwd','AF_rep_U_rev','AF_att_U_fwd','AF_att_U_rev']
for i in range(len(names)):
    name=names[i]
    X=XX[i]

    phase_diagram=PhaseDiagram(model=model)
    phase_diagram.set_directory_name(name)
    phase_diagram.set_indep_vars(X)
    phase_diagram.set_dep_vars_func(func_of_dep_vars)
    phase_diagram.set_plotting_func(plotting_func)
    phase_diagram.set_png()
    if EXECUTE:
        phase_diagram.execute()
    else:
        phase_diagram.plot()

names=['Delta_rep_U_fwd','Delta_rep_U_rev','Delta_att_U_fwd','Delta_att_U_rev','Mag_rep_U_fwd','Mag_rep_U_rev','Mag_att_U_fwd','Mag_att_U_rev','AF_rep_U_fwd','AF_rep_U_rev','AF_att_U_fwd','AF_att_U_rev']
minimised_data = PlotMinimisedData()
minimised_data.set_data_directories(names)
minimised_data.set_directory_name('minimised_data')
minimised_data.set_plotting_func(plotting_func)
minimised_data.minimised_data()
minimised_data.set_png()
minimised_data.plot()
