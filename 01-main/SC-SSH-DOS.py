from lib import *
Uvs=[0, 2,   2.5, 3.8, 1.2, 1.2]
Uws=[0, 1.1, 1.1, 1.1, 0.4, 2.1]
V=20.28
w=1
vv=np.arange(0,4,0.02)
xmin,xmax,ymin,ymax=min(vv),max(vv),-2,2
le,lv=401,len(vv)
mat=np.zeros([4,6,lv,le])

title=[rf'$\overline{{\hat{{\mathcal{{G}}}}(\omega)}}$',r'$\overline{\langle|\hat{M}|\rangle}$',r'$\overline{|\text{Staggered density}|}$',r'$\text{IPR}(|\text{Staggered density}|)$']
TITLES=[r'Mean value of the local density of states. The density of states are plotted as a function of the applied bias energy $\omega$ and intracell, interorbital hopping parameter $v/w$. ', r'Mean absolute value of magnetism. We make the observation that the triplet and Coulomb interactions have the effect of shifting the topological states away from zero-energy, but, importantly, do not shift the transition point away from $v=w$.', r'Mean staggered density. ', r'Inverse participation ratio of the staggered density, as a measure of correlation distance of the charge density wave. ']
save_title=[r'Spin-triplet-SSH-LDOS',r'Spin-triplet-SSH-magnetism',r'Spin-triplet-SSH-staggered',r'Spin-triplet-SSH-IPR']
cmap=['Blues','Reds','Greens','RdPu']

for kk in range(6):
    U_v=Uvs[kk]
    U_w=Uws[kk]
    INITIAL_FIELD='unitary'
    for i,v in enumerate(vv):
        A=Atom([0,0],'A')
        B=Atom([0.25,1],'B')
        A.add_orbital('s')
        B.add_orbital('s')
        lattice_vectors=[[1,0],[0,1]]
        bdg=BogoliubovdeGennes(lattice_vectors,'SSH')
        bdg.add_atom(A)
        bdg.add_atom(B)
        bdg.n_spins=2
        mu=0

        phi=5.32
        zeta_v=3.52
        zeta_w=2.23
        Delta_v=3.92
        Delta_w=2.53

        n_cells=43

        bdg.cut_piece(n_cells, [0])
        bdg.set_onsite(-mu,orbital='s')
        bdg.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B', label='v')
        bdg.set_hopping(-w,hop_vector=[1,0],atom_i='B',atom_f='A', label='w')
        bdg.add_impurities(V,[0,0],label='V')
        #########################################################
        temperature=0
        absolute_convergence_factor=0.00001
        friction=0.7
        max_iterations=1000

        bdg.set_temperature(temperature)

        bdg.set_hubbard_u(-U_v*np.eye(2), atom_i='A', atom_f='B', hop_vector=[0,0])
        hubbard_R=U_w*np.eye(2)
        bdg.set_hubbard_u(hubbard_R, atom_i='B', atom_f='A', hop_vector=[1,0])

        spin_sector=np.eye(2)
        spin_sector[[1,1]]=0

        def set_mean_fields(state):

            if state=='INT':
                bdg.set_hartree([1.25*phi,0.75*phi])

                bdg.set_fock(zeta_v*spin_sector, atom_i='A', atom_f='B', hop_vector=[0,0])
                bdg.set_fock(zeta_w*np.eye(2), atom_i='B', atom_f='A', hop_vector=[1,0])

                bdg.set_gorkov(Delta_v*spin_sector, atom_i='A', atom_f='B', hop_vector=[0,0])
                bdg.set_gorkov(Delta_w*np.eye(2), atom_i='B', atom_f='A', hop_vector=[1,0])

            elif state=='unitary':
                bdg.set_hartree([1.25*phi,0.75*phi])

                bdg.set_fock(zeta_v*np.eye(2), atom_i='A', atom_f='B', hop_vector=[0,0])
                bdg.set_fock(zeta_w*np.eye(2), atom_i='B', atom_f='A', hop_vector=[1,0])

                bdg.set_gorkov(Delta_v*np.eye(2), atom_i='A', atom_f='B', hop_vector=[0,0])
                bdg.set_gorkov(Delta_w*np.eye(2), atom_i='B', atom_f='A', hop_vector=[1,0])

            else:
                bdg._hartree=data._hartree
                bdg._fock=data._fock
                bdg._gorkov=data._gorkov
                bdg._hubbard_indices=data._hubbard_indices
                bdg._anomalous_indices=data._anomalous_indices
                bdg.U_entries=data.U_entries

        set_mean_fields(state=INITIAL_FIELD)

        bdg.self_consistent_calculation(friction=friction, max_iterations=max_iterations, absolute_convergence_factor=absolute_convergence_factor)

        energy_interval=np.linspace(ymin,ymax,le)
        resolution=0.1
        bdg.calculate_greens_function(energy_interval,resolution)

        data=bdg
        INITIAL_FIELD='last'

        for j,energy in enumerate(energy_interval):
            ldos=bdg.local_density_of_states(energy=energy, atom=None, orbital=None)
            dos=np.mean(ldos)
            mat[0,kk,i,j]=dos
            stag_den=bdg.mean_abs_staggered_density(atom_i='A', atom_f='B', energy=energy)
            mat[1,kk,i,j]=stag_den
            stag_ipr=bdg.IPR_abs_staggered_density(atom_i='A', atom_f='B', energy=energy)
            mat[2,kk,i,j]=stag_ipr
            mag=bdg.mean_abs_magnetism(energy=energy, atom=None, orbital=None)
            mat[3,kk,i,j]=mag
    ##################################################


for k in range(4):
    fig, axs = plt.subplots(2, 2, sharex='all', sharey='all')
    vmin=np.min(mat[k])
    vmax=np.max(mat[k])
    for ix in range(2):
        for iy in range(3):
            kk=ix*3+iy
            im=axs[ix,iy].imshow(mat[k,kk].T,extent=[xmin,xmax,ymin,ymax], origin='lower',cmap=cmap[k],vmin=vmin,vmax=vmax)
            U_v=Uvs[ix]
            U_w=Uws[iy]
            txt=rf'$U_v={U_v}, U_w={U_w}$'
            text_box = AnchoredText(txt, frameon=True, loc=4, pad=0.5)
            plt.setp(text_box.patch, facecolor='white', alpha=0.5)
            axs[ix,iy].add_artist(text_box)
    fig.set_size_inches(w=LATEX_WIDTH, h=1.6*LATEX_WIDTH) 
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    xlabel=r'$v/w$'
    ylabel=r'$\omega$'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title[k])
    fig.subplots_adjust(right=0.85)
    cax = fig.add_axes([0.88, 0.15, 0.04, 0.7])
    cbar=plt.colorbar(im,cax=cax)
    TITLE=save_title[k]
    TITLE=os.path.join(FIG,TITLE)
    plt.savefig(TITLE+'.pdf', bbox_inches = "tight")
    text=TITLES[k]+rf'Phase diagram of the SSH model in the superconducting state with $\mu/w={mu}$. A large impurity coupling $V={V}$ breaks the periodic boundary conditions.'
    with open(TITLE+'.txt', 'w') as f:
        f.write(text)


