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
    
    if spin_density_wave:
        bdg.set_hartree(rho+rho_shift,atom='Cu',spin='up')
        bdg.set_hartree(rho-rho_shift,atom='Cu',spin='dn')
        bdg.set_hartree(rho-rho_shift,atom='O_x',spin='up')
        bdg.set_hartree(rho+rho_shift,atom='O_x',spin='dn')
        bdg.set_hartree(rho-rho_shift,atom='O_y',spin='up')
        bdg.set_hartree(rho+rho_shift,atom='O_y',spin='dn')

    elif charge_density_wave:
        bdg.set_hartree(rho+rho_shift,atom='Cu')
        bdg.set_hartree(rho-rho_shift,atom='O_x')
        bdg.set_hartree(rho-rho_shift,atom='O_y')
        bdg.set_hartree_antiferromagnetic(rho,atom='O_y')

    else:
        bdg.set_hartree(rho+rho_shift,spin='up')
        bdg.set_hartree(rho-rho_shift,spin='dn')

    bdg.set_fock(phi*1j*Pauli_y)

    bdg.set_gorkov((chi+chi_shift)*1j*Pauli_y,atom_i='Cu',atom_f='Cu')
    bdg.set_gorkov((chi-chi_shift)*1j*Pauli_y,atom_i='O_x',atom_f='O_x')
    bdg.set_gorkov((chi-chi_shift)*1j*Pauli_y,atom_i='O_y',atom_f='O_y')
    _print=False
    ##############################
    ### Add to notes and docs! ###
    ###     U>0 repulsive      ###
    ###     U<0 attractive     ###
    ##############################
    bdg.set_hubbard_u(U*Pauli_x)
    # bdg.set_hubbard_u_impurities(-U*Pauli_x, atom='Cu', impurity_locations=impurity_locations)
    bdg.U=U

    bdg.record_hartree(location=[0,0], atom='Cu', spin='up', _print=_print)
    bdg.record_hartree(location=[0,0], atom='Cu', spin='down', _print=_print)
    bdg.record_hartree(location=[0,0], atom='O_x', spin='up', _print=_print)
    bdg.record_hartree(location=[0,0], atom='O_x', spin='down', _print=_print)
    bdg.record_fock(location_i=[0,0], location_f=[0,0],atom_i='Cu', atom_f='Cu', spin_i='up', spin_f='dn',_print=_print)
    bdg.record_fock(location_i=[0,0], location_f=[0,0],atom_i='O_x', atom_f='O_x', spin_i='up', spin_f='dn',_print=_print)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='Cu', atom_f='Cu', spin_i='up', spin_f='dn',_print=_print)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='O_x', atom_f='O_x', spin_i='up', spin_f='dn',_print=_print)

    # bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='Cu', atom_f='Cu', spin_i='dn', spin_f='up',_print=_print)
    # bdg.record_gorkov(location_i=[0,0], location_f=[0,0], atom_i='O_x', atom_f='O_x', spin_i='dn', spin_f='up',_print=_print)
    return bdg

def func_of_dep_vars(bdg):
    energy_interval=np.linspace(-8,8,600)
    resolution=0.1
    bdg.greens_function_xy=GreensFunction(bdg,energy_interval,resolution, k_axes=None)

    return bdg

def plotting_func(bdg):
    greens_function = bdg.greens_function_xy

    # Hartree number
    # ==============
    H_Cu = bdg.hartree(atom='Cu')
    H_O_x = bdg.hartree(atom='O_x')
    H_O_y = bdg.hartree(atom='O_x')

    # Hartree magnetism
    # =================
    M_Cu = bdg.hartree_magnetism(atom='Cu')
    M_O_x = bdg.hartree_magnetism(atom='O_x')
    M_O_y = bdg.hartree_magnetism(atom='O_x')
    
    # Gorkov
    # ======
    G_Cu = bdg.gorkov(atom='Cu',spin_i='up',spin_f='dn')
    G_O_x = bdg.gorkov(atom='O_x',spin_i='up',spin_f='dn')
    G_O_y = bdg.gorkov(atom='O_x',spin_i='up',spin_f='dn')

    # DOS
    # ===
    energy = greens_function.energy_interval
    DOS_Cu = greens_function.spectrum(sites='integrated', atom='Cu',orbital='integrated', spin='integrated', anomalous=False)
    DOS_O_x = greens_function.spectrum(sites='integrated', atom='O_x',orbital='integrated', spin='integrated', anomalous=False)
    DOS_O_y = greens_function.spectrum(sites='integrated', atom='O_y',orbital='integrated', spin='integrated', anomalous=False)

    fig,axs=plt.subplots(4,1)

    # Plotting
    # ========
    linestyles=['solid','dashed','solid']

    # Hartree
    axs[0].plot(H_Cu[0],c='r', linestyle=linestyles[0])
    axs[0].plot(H_O_x[0],c='g',linestyle=linestyles[1])
    #axs[0,0].imshow(H_Cu)
    #axs[0,1].imshow(H_O_x)
    #axs[0,2].imshow(H_O_y)
    # Hartree magnetism
    axs[1].plot(M_Cu[0],c='r', linestyle=linestyles[0])
    axs[1].plot(M_O_x[0],c='g',linestyle=linestyles[1])
    #axs[1,0].imshow(M_Cu)
    #axs[1,1].imshow(M_O_x)
    #axs[1,2].imshow(M_O_y)
    # Gorkov
    axs[2].plot(G_Cu[0],c='r', linestyle=linestyles[0])
    axs[2].plot(G_O_x[0],c='g',linestyle=linestyles[1])
    #axs[2,0].imshow(G_Cu)
    #axs[2,1].imshow(G_O_x)
    #axs[2,2].imshow(G_O_y)
    # Dos
    # gs = axs[-1, 1].get_gridspec()
    # # remove the underlying axes
    # for ax in axs[-1, :]:
    #     ax.remove()
    # axbig = fig.add_subplot(gs[-1,:])
    axs[-1].plot(energy,DOS_Cu,label=r'$\text{Cu}^{d_{x^2-y^2}}$',c='r',linestyle=linestyles[0])
    axs[-1].plot(energy,DOS_O_x,label=r'$\text{O}^{p_x}$',c='g',linestyle=linestyles[1],zorder=2)
    axs[-1].plot(energy,DOS_O_y,label=r'$\text{O}^{p_y}$',c='b',linestyle=linestyles[2],zorder=1)
    axs[-1].set_xlabel(r'$\omega$')
    axs[-1].set_ylabel(r'$\text{DOS}(\omega)$')
    axs[-1].legend()


    fig.suptitle(rf'$U={bdg.U:0.03f}$')

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    plt.tight_layout()
    return fig

def phase_diagram_extract_func(bdg):

    greens_function = bdg.greens_function_xy

    N1=bdg.hartree(atom='Cu')[0,0]
    N2=bdg.hartree(atom='O_x')[0,0]
    N3=bdg.hartree(atom='O_y')[0,0]
    M1=np.abs(bdg.hartree(atom='Cu',spin='up')-bdg.hartree(atom='Cu',spin='dn'))[0,0]
    M2=np.abs(bdg.hartree(atom='O_x',spin='up')-bdg.hartree(atom='O_x',spin='dn'))[0,0]
    M3=np.abs(bdg.hartree(atom='O_y',spin='up')-bdg.hartree(atom='O_y',spin='dn'))[0,0]
    Phi1=bdg.fock(atom_i='Cu',atom_f='Cu',spin_i='up',spin_f='dn')[0,0]
    Phi2=bdg.fock(atom_i='O_x',atom_f='O_x',spin_i='up',spin_f='dn')[0,0]
    Phi3=bdg.fock(atom_i='O_y',atom_f='O_y',spin_i='up',spin_f='dn')[0,0]
    Delta1=bdg.gorkov(atom_i='Cu',atom_f='Cu',spin_i='up',spin_f='dn')[0,0]
    Delta2=bdg.gorkov(atom_i='O_x',atom_f='O_x',spin_i='up',spin_f='dn')[0,0]
    Delta3=bdg.gorkov(atom_i='O_y',atom_f='O_y',spin_i='up',spin_f='dn')[0,0]
    N=N1+N2
    CDW=np.abs(N1-N2)
    M=np.abs(M1+M2)
    AF=np.abs(M1-M2)
    Phi=Phi1+Phi2
    PhiDW=np.abs(Phi1-Phi2)
    Delta=np.abs(Delta1+Delta2)
    DeltaDW=np.abs(Delta1-Delta2)

    FE=bdg.free_energy[-1]
    y=[N,M,AF,CDW,Phi,PhiDW,Delta,DeltaDW,FE]
    return y

def phase_diagram_plot_func(x,y):

    fig, [ax1, ax2] = plt.subplots(2,1,sharex='col')

    labels = [
            r'$\langle\hat{\mathcal{N}}\rangle$',
            r'$|\langle\hat{\mathcal{M}}\rangle|$',
            r'$|\text{AF}|$',
            r'$|\text{CDW}|$',
            r'$\Phi$',
            r'$|\Phi\text{DW}|$',
            r'$|\Delta|$',
            r'$|\Delta\text{DW}|$']
    
    markers=['.','x','s','+','^','>','v','o']

    colors=['b','g','r','c','m','y','k','orange']

    ax2.set_ylabel(r'Free energy')
    fig.suptitle(r'Phase diagram')
    fig.supxlabel(r'$U$')
    
    ax1.set_ylabel('Field amplitude')

    for i in range(len(labels)):

        ax1.plot(x,y[i],label=labels[i],marker=markers[i],c=colors[i],markersize=s)

    ax2.plot(x,y[-1],marker='o',c='r',markersize=s)

    ax1.legend()
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    plt.tight_layout()
    return fig,[ax1,ax2]


def Iterations(bdg):
    phiUpCu=bdg._hartree_iterations[0]
    phiDnCu=bdg._hartree_iterations[1]
    phiUpOx=bdg._hartree_iterations[2]
    phiDnOx=bdg._hartree_iterations[3]
    PhiCu=bdg._fock_iterations[0]
    PhiOx=bdg._fock_iterations[1]
    DeltaCu=bdg._gorkov_iterations[0]
    DeltaOx=bdg._gorkov_iterations[1]

    N=phiUpCu+phiDnCu+phiUpOx+phiDnOx
    M=phiUpCu-phiDnCu+phiUpOx-phiDnOx
    AF=(phiUpCu-phiDnCu)-(phiUpOx-phiDnOx)
    CDW=phiUpCu+phiDnCu-(phiUpOx+phiDnOx)
    PHI=PhiCu+PhiOx
    PHIDW=PhiCu-PhiOx
    DELTA=DeltaCu+DeltaOx
    DELTADW=DeltaCu-DeltaOx

    y=[N,M,AF,CDW,PHI,PHIDW,DELTA,DELTADW]

    markers=['.','x','s','+','^','>','v','o']

    colors=['b','g','r','c','m','y','k','orange']

    fig, [ax1, ax2] = plt.subplots(2,1,sharex='col')
    
    s=1

    labels = [
            r'$\langle\hat{\mathcal{N}}\rangle$',
            r'$\langle\hat{\mathcal{M}}\rangle$',
            r'$\text{AF}$',
            r'$\text{CDW}$',
            r'$\Phi$',
            r'$\Phi\text{DW}$',
            r'$\Delta$',
            r'$\Delta\text{DW}$']
    
    for i in range(len(labels)):

        ax1.plot(y[i],label=labels[i],marker=markers[i],c=colors[i],markersize=s)

    ax1.legend()
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    plt.tight_layout()

    ax2.plot(bdg.free_energy,c='r',marker=markers[4],markersize=s)
    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    ax1.set_ylabel('Field amplitude')
    ax2.set_ylabel(r'Free energy')
    fig.suptitle(rf'Self-consistent mean-fields U={bdg.U:.03f}')
    fig.supxlabel(r'Iterations')
    
    plt.tight_layout()

    return fig,[ax1,ax2]


####################################################

n_cells=11

mu=0
s=0
t=1

V=0.00001
L=1

INITIAL=False
EXECUTE=True
PLOT=False

pos_1=np.linspace(0,5,20)
neg_1=np.linspace(-5,0,20)
pos_2=np.linspace(5,10,20)
neg_2=np.linspace(-10,-5,20)
XX=[pos_1[::-1],neg_1,pos_2[::-1],neg_2]

# rho,rho_shift,phi,chi,chi_shift=-2.4,0,0,9.5,0
# charge_density_wave=False
# spin_density_wave=False
# names=['Delta_rep_U_1','Delta_att_U_1','Delta_rep_U_2','Delta_att_U_2']
# for i in range(len(names)):
#     if i==1:
#         rho=10
#     name=names[i]
#     X=XX[i]

#     phase_diagram=PhaseDiagram(model=model)
#     phase_diagram.set_directory_name(name)
#     phase_diagram.set_indep_vars(X)
#     phase_diagram.set_dep_vars_func(func_of_dep_vars)
#     phase_diagram.set_plotting_func(plotting_func)
#     phase_diagram.set_png()
#     phase_diagram.set_plot_initial_iterations(Iterations)
#     if INITIAL:
#         phase_diagram.plot_initial_only()
#         phase_diagram.execute()
#     elif EXECUTE:
#         phase_diagram.execute()
#     elif PLOT:
#         phase_diagram.plot()
#     phase_diagram.set_phase_diagram_extract_func(phase_diagram_extract_func)
#     phase_diagram.set_phase_diagram_plot_func(phase_diagram_plot_func)
#     phase_diagram.plot_phase_diagram()

# rho,rho_shift,phi,chi,chi_shift=-2.6,0,0,4.5,2.2
# charge_density_wave=False
# spin_density_wave=False
# names=['DeltaDW_rep_U_1','DeltaDW_att_U_1','DeltaDW_rep_U_2','DeltaDW_att_U_2']
# for i in range(len(names)):
#     if i==1:
#         rho=10
#     name=names[i]
#     X=XX[i]

#     phase_diagram=PhaseDiagram(model=model)
#     phase_diagram.set_directory_name(name)
#     phase_diagram.set_indep_vars(X)
#     phase_diagram.set_dep_vars_func(func_of_dep_vars)
#     phase_diagram.set_plotting_func(plotting_func)
#     phase_diagram.set_png()
#     phase_diagram.set_plot_initial_iterations(Iterations)
#     if INITIAL:
#         phase_diagram.plot_initial_only()
#         phase_diagram.execute()
#     elif EXECUTE:
#         phase_diagram.execute()
#     elif PLOT:
#         phase_diagram.plot()
#     phase_diagram.set_phase_diagram_extract_func(phase_diagram_extract_func)
#     phase_diagram.set_phase_diagram_plot_func(phase_diagram_plot_func)
#     phase_diagram.plot_phase_diagram()

# rho,rho_shift,phi,chi,chi_shift=0,10,0,0,0
# charge_density_wave=False
# spin_density_wave=False
# names=['Mag_rep_U_1','Mag_att_U_1','Mag_rep_U_2','Mag_att_U_2']
# for i in range(len(names)):
#     name=names[i]
#     X=XX[i]

#     phase_diagram=PhaseDiagram(model=model)
#     phase_diagram.set_directory_name(name)
#     phase_diagram.set_indep_vars(X)
#     phase_diagram.set_dep_vars_func(func_of_dep_vars)
#     phase_diagram.set_plotting_func(plotting_func)
#     phase_diagram.set_png()
#     phase_diagram.set_plot_initial_iterations(Iterations)
#     if INITIAL:
#         phase_diagram.plot_initial_only()
#         phase_diagram.execute()
#     elif EXECUTE:
#         phase_diagram.execute()
#     elif PLOT:
#         phase_diagram.plot()
#     phase_diagram.set_phase_diagram_extract_func(phase_diagram_extract_func)
#     phase_diagram.set_phase_diagram_plot_func(phase_diagram_plot_func)
#     phase_diagram.plot_phase_diagram()

# charge_density_wave=False
# spin_density_wave=True
# rho,rho_shift,phi,chi,chi_shift=0,4,0,0,0
# names=['AF_rep_U_1','AF_att_U_1','AF_rep_U_2','AF_att_U_2']
# for i in range(len(names)):
#     name=names[i]
#     X=XX[i]

#     phase_diagram=PhaseDiagram(model=model)
#     phase_diagram.set_directory_name(name)
#     phase_diagram.set_indep_vars(X)
#     phase_diagram.set_dep_vars_func(func_of_dep_vars)
#     phase_diagram.set_plotting_func(plotting_func)
#     phase_diagram.set_png()
#     phase_diagram.set_plot_initial_iterations(Iterations)
#     if INITIAL:
#         phase_diagram.plot_initial_only()
#         phase_diagram.execute()
#     elif EXECUTE:
#         phase_diagram.execute()
#     elif PLOT:
#         phase_diagram.plot()
#     phase_diagram.set_phase_diagram_extract_func(phase_diagram_extract_func)
#     phase_diagram.set_phase_diagram_plot_func(phase_diagram_plot_func)
#     phase_diagram.plot_phase_diagram()

charge_density_wave=True
spin_density_wave=False
rho,rho_shift,phi,chi,chi_shift=0,10,0,0,0
names=['CDW_rep_U_1','CDW_att_U_1','CDW_rep_U_2','CDW_att_U_2']
for i in range(len(names)):
    name=names[i]
    X=XX[i]

    phase_diagram=PhaseDiagram(model=model)
    phase_diagram.set_directory_name(name)
    phase_diagram.set_indep_vars(X)
    phase_diagram.set_dep_vars_func(func_of_dep_vars)
    phase_diagram.set_plotting_func(plotting_func)
    phase_diagram.set_png()
    phase_diagram.set_plot_initial_iterations(Iterations)
    if INITIAL:
        phase_diagram.plot_initial_only()
        phase_diagram.execute()
    elif EXECUTE:
        phase_diagram.execute()
    elif PLOT:
        phase_diagram.plot()
    phase_diagram.set_phase_diagram_extract_func(phase_diagram_extract_func)
    phase_diagram.set_phase_diagram_plot_func(phase_diagram_plot_func)
    phase_diagram.plot_phase_diagram()

names=['DeltaDW_rep_U_1','DeltaDW_att_U_1','Delta_rep_U_1','Delta_att_U_1','Mag_rep_U_1','Mag_att_U_1','AF_rep_U_1','AF_att_U_1','CDW_rep_U_1','CDW_att_U_1','DeltaDW_rep_U_2','DeltaDW_att_U_2','Delta_rep_U_2','Delta_att_U_2','Mag_rep_U_2','Mag_att_U_2','AF_rep_U_2','AF_att_U_2','CDW_rep_U_2','CDW_att_U_2']
minimised_data = PlotMinimisedData()
minimised_data.set_data_directories(names)
minimised_data.set_directory_name('minimised_data')
minimised_data.set_plotting_func(plotting_func)
minimised_data.minimised_data()

minimised_data.set_png()

if PLOT:
    minimised_data.plot()

minimised_data.set_phase_diagram_extract_func(phase_diagram_extract_func)
minimised_data.set_phase_diagram_plot_func(phase_diagram_plot_func)
minimised_data.plot_phase_diagram()
