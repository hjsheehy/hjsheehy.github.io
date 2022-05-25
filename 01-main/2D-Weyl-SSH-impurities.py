from lib import *

filename=sys.argv[0].split('.')[0]
DATA=os.path.join(DATA,filename+'.npz')
FIG=os.path.join(FIG,filename)
for directory in [FIG]:
    if not os.path.exists(directory):
        os.makedirs(directory)

def main(V,length,glue_edgs,td,v):
    A=Atom([0,0],'A')
    B=Atom([0.5,0],'B')
    A.add_orbital('s')
    B.add_orbital('s')
    lattice_vectors=[[1,0],[0,1]]
    tb=TightBinding(lattice_vectors,'2D-Weyl-SSH')
    tb.add_atom(A)
    tb.add_atom(B)
    tb.n_spins=1
    tb.cut(n_cells, axes=0, glue_edgs=glue_edgs)
    tb.cut(n_cells, axes=1, glue_edgs=True)
    tb.set_onsite(-mu+s,atom='A')
    tb.set_onsite(-mu-s,atom='B')

    tb.set_hopping(-v,hop_vector=[0,0],atom_i='A',atom_f='B',label='$v$')
    tb.set_hopping(-w,hop_vector=[-1,0],atom_i='A',atom_f='B',label='$w$')
    tb.set_hopping(-td,hop_vector=[0,1],atom_i='B',atom_f='A',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[0,-1],atom_i='B',atom_f='A',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[1,1],atom_i='B',atom_f='A',label='$t_d$')
    tb.set_hopping(-td,hop_vector=[1,-1],atom_i='B',atom_f='A',label='$t_d$')

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

            with open(DATA+FIGNAME+'.npz', 'wb') as f:
                cPickle.dump([xx,yy,greens_list], f)
        else:
            [xx,yy,greens_list]=np.load(DATA+FIGNAME+'.npz', allow_pickle=True)

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
        plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight")
        
        if z:
            with open(OUTPUT+'.txt', 'w') as f:
                f.write(rf'''Phase diagram of an impurity wall, perpendicular to the SSH dimers, of length $L$ and coupling strength $V$, and with closed boundary conditions. 
The spectrum of the simple Weyl-SSH model in the normal state, in k-space, integrated over the direction along the SSH chains and summed over the atomic sites.
The model consists of an ${n_cells}\times{n_cells}$ lattice, with chemical potential $\mu={mu}$, interdimer $w={w}$, intradimer $v={v}$ and interchain $t_d={td}$ hoppings.'''+r'''
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

    for i in range(len(FIGNAMES)):
        FIGNAME=FIGNAMES[i]
        glue_edgs=glue_edgss[i]

        r=c=4
        vv=[0,0.6,1.1,1.2,1.3,2.4]
        tdd=[0,0.6,1.1,1.2,1.3,2.4]
        
        if CALCULATE:
            greens_list=[]
            for td in tdd:
                green_list=[]
                for v in vv:
                    greens_function_kq = main(V,length,glue_edgs,td,v)
                    green_list.append(greens_function_kq)
                greens_list.append(green_list)

            with open(DATA+FIGNAME+'.npz', 'wb') as f:
                cPickle.dump([tdd,vv,greens_list], f)
        else:
            [tdd,vv,greens_list]=np.load(DATA+FIGNAME+'.npz', allow_pickle=True)

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
        plt.savefig(OUTPUT+'.pdf', bbox_inches = "tight")
        
        if glue_edgs:
            with open(OUTPUT+'.txt', 'w') as f:
                f.write(rf'''the spectrum of the simple $\text{{nic}}_2$ model in the normal state, in k-space, integrated over the direction along the ssh chains and summed over the atomic sites.
the model consists of an ${n_cells}\times{n_cells}$ lattice, with closed boundary conditions.
the chemical potential is $\mu={mu}$, the interdimer hopping is fixed $w={w}$, and the intradimer $v$ and interchain $t_d$ hoppings are varied.'''+r'''
compare with the model with open boundary conditions fig \ref{fig:phase_diagram_open}.
''')
        else:
            with open(OUTPUT+'.txt', 'w') as f:
                f.write(r'''the model and quantities depicted as in fig \ref{fig:phase_diagram_closed}, but with the boundary conditions closed.
we observe the diagonal hopping $t_d$ has the effect of {\it squeezing} the topological modes into a narrower interval $\mathcal{i}\subset[-\pi,\pi]$. 
at the conventional ssh topological transition point $v=w$, the narrow interval transitions to its complement interval $\mathcal{i}^c$.
notice the universality with the $\text{nic}_2$ model fig \ref{fig:weyl_phase_diagram_open}, \ref{fig:weyl_phase_diagram_closed}.
''')

#############################################################################
################################# Main ######################################
#############################################################################
mu=0.0
s=0.0
td=0.9
v=0.6
w=1.2
n_cells=43

# k_space_phase_diagram(CALCULATE=False)
w_v_phase_diagram(CALCULATE=True)
