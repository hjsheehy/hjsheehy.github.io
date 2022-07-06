from PhaseDiagram import *

def model(U):

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
    
    if antiferromagnetic:
        bdg.set_hartree_antiferromagnetic(rho_shift)
        bdg.set_hartree(rho,spin='up')
        bdg.set_hartree(rho,spin='dn')
    else:
        bdg.set_hartree(rho+rho_shift,spin='up')
        bdg.set_hartree(rho-rho_shift,spin='dn')
        bdg.set_hartree_antiferromagnetic(rho_shift)
    bdg.set_fock(phi*1j*Pauli_y)
    bdg.set_gorkov(chi*1j*Pauli_y)
    
    _print=False
    ##############################
    ### Add to notes and docs! ###
    ###     U>0 repulsive      ###
    ###     U<0 attractive     ###
    ##############################
    bdg.set_hubbard_u(U*Pauli_x)
    bdg.U=U

    bdg.record_hartree(location=[0,0], spin='up', _print=_print)
    bdg.record_hartree(location=[0,0], spin='dn', _print=_print)
    bdg.record_hartree(location=[1,0], spin='up', _print=_print)
    bdg.record_hartree(location=[1,0], spin='dn', _print=_print)
    bdg.record_fock(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='dn',_print=_print)
    bdg.record_gorkov(location_i=[0,0], location_f=[0,0], spin_i='up', spin_f='dn',_print=_print)

    return bdg

def func_of_dep_vars(bdg):
    energy_interval=np.linspace(-8,8,600)
    resolution=0.1
    bdg.greens_function_xy=GreensFunction(bdg,energy_interval,resolution, k_axes=None)

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
    greens_function = bdg.greens_function_xy

    x = greens_function.energy_interval
    y = greens_function.spectrum(sites='integrated', atom='integrated',orbital='integrated', spin='integrated', anomalous=False)
    ax.plot(y,x)
    ax.set_xlabel(r'$\text{DOS}(\omega)$')
    ax.set_ylabel(r'$\omega$')
    N=np.mean(bdg.hartree())
    M1=np.array([(bdg.hartree(spin='up')-bdg.hartree(spin='dn'))[i,i%2::2] for i in range(n_cells)])
    M2=np.array([(bdg.hartree(spin='up')-bdg.hartree(spin='dn'))[i,(i+1)%2::2] for i in range(n_cells)])
    M=M1+M2
    AF=M1-M2
    M=float(np.abs(st.mode(M,axis=None))[0])
    AF=float(np.abs(st.mode(AF,axis=None))[0])
    Phi=np.abs(np.mean(bdg.fock(spin_i='up',spin_f='dn')))
    Delta=np.abs(np.mean(bdg.gorkov(spin_i='up',spin_f='dn')))
    ax.set_title(rf'''$U={bdg.U:0.03f}$
$N={N:0.03f}, M={M:0.03f}, AF={AF:0.03f}, \Delta={Delta:0.03f}$''')

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    plt.tight_layout()
    return fig

def phase_diagram_extract_func(bdg):

    greens_function = bdg.greens_function_xy

    N=np.mean(bdg.hartree())
    M1=np.array([(bdg.hartree(spin='up')-bdg.hartree(spin='dn'))[i,i%2::2] for i in range(n_cells)])
    M2=np.array([(bdg.hartree(spin='up')-bdg.hartree(spin='dn'))[i,(i+1)%2::2] for i in range(n_cells)])
    M=M1+M2
    AF=M1-M2
    # M=np.abs(np.mean(M))
    # AF=np.abs(np.mean(AF))
    M=float(np.abs(st.mode(M,axis=None))[0])
    AF=float(np.abs(st.mode(AF,axis=None))[0])
    Phi=np.abs(np.mean(bdg.fock(spin_i='up',spin_f='dn')))
    Delta=np.abs(np.mean(bdg.gorkov(spin_i='up',spin_f='dn')))
    FE=bdg.free_energy[-1]
    y=[N,M,AF,Phi,Delta,FE]
    return y

def phase_diagram_plot_func(x,y):

    fig, [ax1, ax2] = plt.subplots(2,1,sharex='col')

    labels = [
            r'$\langle\hat{\mathcal{N}}\rangle$',
            r'$|\langle\hat{\mathcal{M}}\rangle|$',
            r'$|\text{AF}|$',
            r'$\Phi$',
            r'$|\Delta|$']
    
    markers=['.','x','s','^','v']

    colors=['b','g','r','m','k']

    ax2.set_ylabel(r'Free energy')
    fig.suptitle(r'Phase diagram')
    fig.supxlabel(r'$U$')
    
    ax1.set_ylabel('Field amplitude')

    s=1

    for i in range(len(labels)):

        ax1.plot(x,y[i],label=labels[i],marker=markers[i],c=colors[i],markersize=s)

    ax2.plot(x,y[-1],marker='o',c='r',markersize=s)

    ax1.legend()

    # ax3 = ax1.twinx()
    # color = 'r'
    # ax3.set_ylabel('Antiferromagnetism', color=color)  # we already handled the x-label with ax1
    # ax3.tick_params(axis='y', labelcolor=color)

    # i=2
    # ax3.plot(x,y[i],label=labels[i],marker=markers[i],c=colors[i],markersize=s)

    fig.set_size_inches(w=LATEX_WIDTH, h=LATEX_WIDTH) 
    plt.tight_layout()
    return fig,[ax1,ax2]

def Iterations(bdg):
    phiUp0=bdg._hartree_iterations[0]
    phiDn0=bdg._hartree_iterations[1]
    phiUp1=bdg._hartree_iterations[2]
    phiDn1=bdg._hartree_iterations[3]
    Phi=bdg._fock_iterations[0]
    Delta=bdg._gorkov_iterations[0]

    N=(phiUp0+phiDn0+phiUp1+phiDn1)/2
    M0=(phiUp0-phiDn0)/2
    M1=(phiUp1-phiDn1)/2
    M=(M0+M1)/2
    AF=(M0-M1)/2
    PHI=Phi
    DELTA=Delta

    y=[N,M,AF,PHI,DELTA]

    markers=['.','x','s','^','v']

    colors=['b','g','r','m','k']

    fig, [ax1, ax2] = plt.subplots(2,1,sharex='col')
    
    s=1

    labels = [
            r'$\langle\hat{\mathcal{N}}\rangle$',
            r'$\langle\hat{\mathcal{M}}\rangle$',
            r'$\text{AF}$',
            r'$\Phi$',
            r'$\Delta$']
    
    
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

n_cells=8
mu=0
t=1

V=0
impurity_locations=[[0,0]]

antiferromagnetic=False

pos=np.linspace(-10,10,81)
neg=pos[:61]
pos=pos[20:]
pos2=pos[61::-1]
pos1=pos[38:]
neg1=neg[:41]
neg2=neg[30::-1]
XX=[pos1,pos2,neg1,neg2]

# rho,rho_shift,phi,chi=5,0,0,4
# names=['Delta_att_U_fwd','Delta_att_U_rev']
# for i in range(len(names)):
#     name=names[i]
#     i+=2
#     X=XX[i]

#     phase_diagram=PhaseDiagram(model=model)
#     phase_diagram.set_directory_name(name)
#     phase_diagram.set_indep_vars(X)
#     phase_diagram.set_dep_vars_func(func_of_dep_vars)
#     phase_diagram.set_plotting_func(plotting_func)
#     phase_diagram.set_phase_diagram_plot_func(phase_diagram_plot_func)
#     phase_diagram.set_plot_initial_iterations(Iterations)
#     # phase_diagram.execute()
#     # phase_diagram.set_png()
#     phase_diagram.plot()
#     phase_diagram.set_phase_diagram_extract_func(phase_diagram_extract_func)
#     phase_diagram.plot_phase_diagram()

# rho,rho_shift,phi,chi=10,10,0,0
# names=['Mag_rep_U_fwd','Mag_rep_U_rev']
# for i in range(len(names)):
#     i=1
#     name=names[i]
#     X=XX[i]

#     phase_diagram=PhaseDiagram(model=model)
#     phase_diagram.set_directory_name(name)
#     phase_diagram.set_indep_vars(X)
#     phase_diagram.set_dep_vars_func(func_of_dep_vars)
#     phase_diagram.set_plotting_func(plotting_func)
#     phase_diagram.set_phase_diagram_plot_func(phase_diagram_plot_func)
#     phase_diagram.set_plot_initial_iterations(Iterations)
#     phase_diagram.execute()
#     # phase_diagram.set_png()
#     phase_diagram.plot()
#     phase_diagram.set_phase_diagram_extract_func(phase_diagram_extract_func)
#     phase_diagram.plot_phase_diagram()
#     exit()


rho,rho_shift,phi,chi=2,2,0,0
names=['AF_rep_U_fwd','AF_rep_U_rev']
antiferromagnetic=True

rho=0
bdg=model(4)
bdg.self_consistent_calculation(friction=0.7, max_iterations=1000, absolute_convergence_factor=0.0001)
energy_interval=np.linspace(-8,8,600)
resolution=0.1
greens_function_xy=GreensFunction(bdg,energy_interval,resolution, k_axes=None)
mag=greens_function_xy.magnetism(sites='resolved', atom='integrated', orbital='integrated', energy=0, anomalous=False)
mag=np.fft.fftshift(mag)
im=plt.imshow(mag)
plt.colorbar(im)
plt.show()

N=np.mean(bdg.hartree())
M1=np.array([(bdg.hartree(spin='up')-bdg.hartree(spin='dn'))[i,i%2::2] for i in range(n_cells)])
M2=np.array([(bdg.hartree(spin='up')-bdg.hartree(spin='dn'))[i,(i+1)%2::2] for i in range(n_cells)])
M=M1+M2
AF=M1-M2
# M=np.abs(np.mean(M))
# AF=np.abs(np.mean(AF))
M=float(np.abs(st.mode(M,axis=None))[0])
AF=float(np.abs(st.mode(AF,axis=None))[0])
Phi=np.abs(np.mean(bdg.fock(spin_i='up',spin_f='dn')))
Delta=np.abs(np.mean(bdg.gorkov(spin_i='up',spin_f='dn')))
exit()

for i in range(len(names)):
    if i==1:
        rho=0
        rho_shift=10
    name=names[i]
    X=XX[i]

    phase_diagram=PhaseDiagram(model=model)
    phase_diagram.set_directory_name(name)
    phase_diagram.set_indep_vars(X)
    phase_diagram.set_dep_vars_func(func_of_dep_vars)
    phase_diagram.set_plotting_func(plotting_func)
    phase_diagram.set_phase_diagram_plot_func(phase_diagram_plot_func)
    phase_diagram.set_plot_initial_iterations(Iterations)
    phase_diagram.execute()
    # phase_diagram.set_png()
    phase_diagram.plot()
    phase_diagram.set_phase_diagram_extract_func(phase_diagram_extract_func)
    phase_diagram.plot_phase_diagram()

exit()
names=['Delta_rep_U_fwd','Delta_rep_U_rev','Delta_att_U_fwd','Delta_att_U_rev','Mag_rep_U_fwd','Mag_rep_U_rev','Mag_att_U_fwd','Mag_att_U_rev','AF_rep_U_fwd','AF_rep_U_rev','AF_att_U_fwd','AF_att_U_rev']
minimised_data = PlotMinimisedData()
minimised_data.set_data_directories(names)
minimised_data.set_directory_name('minimised_data')
minimised_data.set_plotting_func(plotting_func)
minimised_data.minimised_data()
# minimised_data.set_png()
minimised_data.plot()
minimised_data.set_phase_diagram_extract_func(phase_diagram_extract_func)
minimised_data.set_phase_diagram_plot_func(phase_diagram_plot_func)
minimised_data.plot_phase_diagram()
