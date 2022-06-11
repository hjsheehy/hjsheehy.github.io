from lib import *

FILENAME=sys.argv[0].split('.')[0]
DATA=os.path.join(DATA,FILENAME)
FIG=os.path.join(FIG,FILENAME)
for directory in [FIG,DATA]:
    if not os.path.exists(directory):
        os.makedirs(directory)

def main(V,length,glue_edgs,td,v):
    A=Atom([0,0],'A')
    B=Atom([0.5,0],'B')
    A.add_orbital('s')
    B.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    tb=TightBinding(lattice_vectors,'LaNiC2')
    tb.add_atom(A)
    tb.add_atom(B)
    tb.n_spins=1
    tb.cut(nx, axes=0, glue_edgs=glue_edgs)
    tb.cut(ny, axes=1, glue_edgs=True)
    tb.set_onsite(-mu+s,atom='A')
    tb.set_onsite(-mu-s,atom='B')

    tb.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B',label='$v$')
    tb.set_hopping(-w,hop_vector=[-1,0],atom_i='A',atom_f='B',label='$w$')
    tb.set_hopping(-td,hop_vector=[0,1],atom_i='B',atom_f='A',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[0,-1],atom_i='B',atom_f='A',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[-1,1],atom_i='A',atom_f='B',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[-1,-1],atom_i='A',atom_f='B',label='$t_d$')

    locations=[[0,0]]
    for l in range(1,length):
        locations.append([0,-l])
        locations.append([0,l])
    tb.add_impurities(V, locations)

    tb.solve()

    energy_interval=np.linspace(-4,4,601)
    resolution=0.05

    greens_function_kq=GreensFunction(tb,energy_interval,resolution, k_axes=[0,1])

    del tb.eigenvectors
    del tb.eigenvalues

    return greens_function_kq

# Plotting:
def k_space_phase_diagram(CALCULATE):
    
    FIGNAMES=['phase_diagram_closed','phase_diagram_open']
    glue_edgss=[True,False]

    x_label='V'
    xx=[0,0.1,0.5,1.2,1.9,5,100]
    y_label='L'
    yy=[0,1,3,5,7,9,21]

    for i in range(2):
        FIGNAME=FIGNAMES[i]
        z=glue_edgss[i]

        r=c=4
        
        if CALCULATE:
            greens_list=[]
            for x in xx:
                green_list=[]
                for y in yy:
                    greens_function_kq = main(x,y,z,td,v)
                    green_list.append(greens_function_kq)
                greens_list.append(green_list)

            with open(os.path.join(DATA,FIGNAME+'.npz'), 'wb') as f:
                cPickle.dump([xx,yy,greens_list], f)
        else:
            [xx,yy,greens_list]=np.load(os.path.join(DATA,FIGNAME+'.npz'), allow_pickle=True)

        fig, axs = plt.subplots(len(yy),len(xx))

        bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
        for i in range(len(yy)):
            for j in range(len(xx)):
                k=len(yy)-i-1
                greens_function = greens_list[j][k]
                x,y=xx[j],yy[k]
                axs[i,j] = greens_function.plot_spectrum(axs[i,j], axes=['integrated','resolved'],omega_min=-2,omega_max=2,vmin=0,vmax=20, atom='integrated')
                axs[i,j].set_title('')
                axs[i,j].set_xlabel('')
                axs[i,j].set_ylabel('')
                if i!=len(yy)-1:
                    axs[i,j].set_xticks([])
                if j!=0:
                    axs[i,j].set_yticks([])
                if j==len(xx)-1:
                    axs[i,j].yaxis.set_label_position("right")
                    axs[i,j].set_ylabel(rf'${2*y+1}$')
                    if i==len(yy)-1:
                        axs[i,j].set_ylabel(rf'${y_label}={2*y+1}$')
                if i==0:
                    axs[i,j].xaxis.set_label_position("top")
                    axs[i,j].set_xlabel(rf'${x}$')
                    if j==0:
                        axs[i,j].set_xlabel(rf'${x_label}={x}$')
                        

        fig.suptitle(greens_function.title)
        fig.supxlabel(greens_function.xlabel)
        fig.supylabel(greens_function.ylabel)
        fig.set_size_inches(w=LATEX_WIDTH, h=1.2*LATEX_WIDTH) 
        plt.subplots_adjust(wspace=0.3, hspace=0.25)
        # plt.tight_layout()

        OUTPUT=os.path.join(FIG,FIGNAME)
        plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight",dpi=DPI)
        
        if z:
            with open(OUTPUT+'.txt', 'w') as f:
                f.write(rf'''Phase diagram of an impurity wall, perpendicular to the SSH dimers, of length $L$ and coupling strength $V$, and with closed boundary conditions. 
The spectrum of the simple Weyl-SSH model in the normal state, in k-space, integrated over the direction along the SSH chains and summed over the atomic sites.
The model consists of an ${nx}\times{ny}$ lattice, with chemical potential $\mu={mu}$, interdimer $w={w}$, intradimer $v={v}$ and interchain $t_d={td}$ hoppings.'''+r'''
Compare with the model with open boundary conditions Fig \ref{fig:wall_phase_diagram_open}.
''')
        else:
            with open(OUTPUT+'.txt', 'w') as f:
                f.write(r'''The model and quantities depicted as in Fig \ref{fig:wall_phase_diagram_closed}, but with the boundary conditions open.
We observe the topological mode due to the impurity is independent of the topological modees assoicated with the boundary conditions.
''')

def w_v_phase_diagram(CALCULATE):
    
    FIGNAMES=['w_v_phase_diagram']
    glue_edgss=[True]
    V=100
    length=int(n_cells/2)

    for ii in range(len(FIGNAMES)):
        FIGNAME=FIGNAMES[ii]
        glue_edgs=glue_edgss[ii]

        r=c=4
        vv=[0,0.3,0.6,1.1,1.2,1.3,2.4,5,10,20]
        tdd=[0,0.3,0.6,1.1,1.2,1.3,2.4,5,10,20]
        
        if CALCULATE:
            greens_list=[]
            for td in tdd:
                green_list=[]
                for v in vv:
                    greens_function_kq = main(V,length,glue_edgs,td,v)
                    green_list.append(greens_function_kq)
                greens_list.append(green_list)

            with open(os.path.join(DATA,FIGNAME)+'.npz', 'wb') as f:
                cPickle.dump([tdd,vv,greens_list], f)
        else:
            [tdd,vv,greens_list]=np.load(os.path.join(DATA,FIGNAME+'.npz'), allow_pickle=True)

        fig, axs = plt.subplots(len(vv),len(tdd))

        bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
        for i in range(len(vv)):
            for j in range(len(tdd)):
                k=len(vv)-i-1
                greens_function = greens_list[j][k]
                td,v=tdd[j],vv[k]
                axs[i,j] = greens_function.plot_spectrum(axs[i,j], axes=['integrated','resolved'],omega_min=-2,omega_max=2,vmin=0,vmax=20, atom='integrated')
                axs[i,j].set_title('')
                axs[i,j].set_xlabel('')
                axs[i,j].set_ylabel('')
                if i!=len(vv)-1:
                    axs[i,j].set_xticks([])
                if j!=0:
                    axs[i,j].set_yticks([])
                if j==len(tdd)-1:
                    axs[i,j].yaxis.set_label_position("right")
                    axs[i,j].set_ylabel(rf'${v}$')
                    if i==len(vv)-1:
                        axs[i,j].set_ylabel(rf'$v/w={v}$')
                if i==0:
                    axs[i,j].xaxis.set_label_position("top")
                    axs[i,j].set_xlabel(rf'${td}$')
                    if j==0:
                        axs[i,j].set_xlabel(rf'$t_d={td}$')
                        

        fig.suptitle(greens_function.title)
        fig.supxlabel(greens_function.xlabel)
        fig.supylabel(greens_function.ylabel)
        fig.set_size_inches(w=LATEX_WIDTH, h=1.2*LATEX_WIDTH) 
        plt.subplots_adjust(wspace=0.3, hspace=0.25)
        # plt.tight_layout()

        OUTPUT=os.path.join(FIG,FIGNAME)
        plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight",dpi=DPI)
        
        if glue_edgs:
            with open(OUTPUT+'.txt', 'w') as f:
                f.write(rf'''The spectrum of the simple $\text{{NiC}}_2$ model in the normal state, in k-space, integrated over the direction along the SSH chains and summed over the atomic sites.
The model consists of an ${nx}\times{ny}$ lattice, with closed boundary conditions.
The chemical potential is $\mu={mu}$, the interdimer hopping is fixed $w={w}$, and the intradimer $v$ and interchain $t_d$ hoppings are varied.'''+r'''
Compare with the model with open boundary conditions Figure \ref{fig:phase_diagram_open}.
''')
        else:
            with open(OUTPUT+'.txt', 'w') as f:
                f.write(r'''The model and quantities depicted as in Figure \ref{fig:phase_diagram_closed}, but with the boundary conditions closed.
We observe the diagonal hopping $t_d$ has the effect of creating a boundary in the Brillouin zone, between its centre and its outside, narrowing the centre with increasing $t_d$.
In the conventional SSH model, the topological transition point, $v=w$, is a metallic phase, between a trivial insulator and a Majorana topological insulator.
In this 2D model, the transition point $w=v$ is a topological insulator with a pair of Majorana Dirac cones {\color{red}check this}, which move towards the centre of the Brillouin zone with increasing $t_d$, a second order transition, with the Dirac nodes merging at $t_d\to\infty$.
The inside of the Brillouin zone, between the Dirac cones, exhibits a Fermi arc in the SSH topological phase, $v<w$.
The transition, $v=w$, is a first-order transition; the Fermi arc splits into two, one sharply moving up the spectrum, into the bulk, as the other moves downwards, which is a Majorana topological insulting phase.
The SSH trivial phase, $v>w$, becomes non-trivial in 2D; with the Fermi arc in the opposite side of the Brillouin zone.
Notice the universality with the $\text{NiC}_2$ model Figure \ref{fig:weyl_phase_diagram_open}, \ref{fig:weyl_phase_diagram_closed}.
''')

def V_phase_diagram(CALCULATE):
    
    FIGNAMES=['V_phase_diagram_topological','V_phase_diagram_transition','V_phase_diagram_trivial']
    glue_edgss=[True,True,True]
    vv=[0.6,1.2,1.8]
    VV=[-100,-40,-20,-10,-8.0,-7.5,-7.0,-6.5,-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,
            0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,10,20,40,100]
    VV=VV[::-1]
    length=int(n_cells/2)

    for ii in range(len(FIGNAMES)):
        FIGNAME=FIGNAMES[ii]
        glue_edgs=glue_edgss[ii]
        v=vv[ii]

        r=c=4
        
        if CALCULATE:
            greens_list=[]
            for V in VV:
                greens_function_kq = main(V,length,glue_edgs,td,v)
                greens_list.append(greens_function_kq)
            with open(os.path.join(DATA,FIGNAME+'.npz'), 'wb') as f:
                cPickle.dump([VV,greens_list], f)
        else:
            [VV,greens_list]=np.load(os.path.join(DATA,FIGNAME+'.npz'), allow_pickle=True)
        
        N=int(np.sqrt(len(VV)))
        fig, axs = plt.subplots(N,N)

        bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
        k=0
        for i in range(N):
            for j in range(N):
                greens_function = greens_list[k]
                V=VV[k]
                axs[i,j] = greens_function.plot_spectrum(axs[i,j], axes=['integrated','resolved'],omega_min=-2,omega_max=2,vmin=0,vmax=10, atom='integrated')
                axs[i,j].set_title('')
                axs[i,j].set_xlabel('')
                axs[i,j].set_ylabel('')
                props = dict(boxstyle='round', facecolor='w', alpha=0.5)
                textstr=rf'${V}$'
                if k==0:
                    textstr=rf'$V={V}$'
                # place a text box in upper left in axes coords
                axs[i,j].text(0.5, 0.9, textstr, transform=axs[i,j].transAxes,
                        verticalalignment='top', horizontalalignment='center', bbox=props)

                if i!=N-1:
                    axs[i,j].set_xticks([])
                if j!=0:
                    axs[i,j].set_yticks([])

                k+=1

        fig.suptitle(greens_function.title)
        fig.supxlabel(greens_function.xlabel)
        fig.supylabel(greens_function.ylabel)
        fig.set_size_inches(w=LATEX_WIDTH, h=1.2*LATEX_WIDTH) 
        plt.subplots_adjust(wspace=0.3, hspace=0.25)
        # plt.tight_layout()

        OUTPUT=os.path.join(FIG,FIGNAME)
        # plt.show()
        # exit()
        plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight",dpi=DPI)
        
        if ii==0:
            with open(OUTPUT+'.txt', 'w') as f:
                f.write(rf'''Impurity wall. The spectrum of the Weyl-SSH model in the normal state, in k-space, integrated over the direction along the SSH chains and summed over the atomic sites, ${nx}\times{ny}$ cells, with closed boundary conditions.
    The chemical potential is $\mu={mu}$, interdimer hopping $w={w}$, intradimer $v={v}$ and interchain $t_d={td}$, corresponding to the topological phase of the SSH model.
    The potential barrier is attractive with negative $V<0$ and repulsive with positive $V>0$.'''+r'''
    For small $|V|\approx w$, a mode outside the Brillouin Zone leaves the bulk.
    In particular, zero-energy modes exist at the edge of the Brillouin zone when $|V|=w$.
    As $|V|>w$ increases, the mode moves into the centre of the Brillouin Zone, falling or rising with the sign of $V$, tending to the topological spectrum of the open-boundary model as $|V|\to\infty$.
    ''')
        elif ii==1:
            with open(OUTPUT+'.txt', 'w') as f:
                f.write(r'''Model and quantities depicted in Figure \ref{fig:V_phase_diagram_topological}, except with
'''+rf'''the interdimer hopping $v=w={v}$, the topological transition point of the SSH model.
''')
        elif ii==2:
            with open(OUTPUT+'.txt', 'w') as f:
                f.write(r'''Model and quantities depicted in Figure \ref{fig:V_phase_diagram_topological}, except with
'''+rf'''the interdimer hopping $v={v}>w$, the trivial phase of the SSH model.
''')
#############################################################################
################################# Main ######################################
#############################################################################
mu=0.0
s=0.0
td=0.9
v=0.6
w=1.2
n_cells=nx=ny=43

k_space_phase_diagram(CALCULATE=True)
w_v_phase_diagram(CALCULATE=True)
V_phase_diagram(CALCULATE=True)
