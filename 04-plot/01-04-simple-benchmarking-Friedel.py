from plt_lib import *
########################################################
########################## Data ########################
########################################################
layer=0
omega=0
n_data=len(filenames)
n_list=[]
dis_straight=[]
dos_straight=[]
dis_diagonal=[]
dos_diagonal=[]
dis_imp=[]

for i in range(n_data):
    ################ import ################
    filename = filenames[i]
    globals().update(conf_file(filename))

    data = np.load(filename, allow_pickle=True)
    dos, exec_time, mem = data['dos'], data['exec_time'], data['mem']
    ldos = LDOS(dos, omegas, omega, trace_over=True, layer=layer)
    
    ###### ldos from impurity to edge ######
    model=Model(dimensions, n_spins, basis, orbitals, hop_links, impurities)
    centre=model.centre
    y=0
    n_list.append(n_x)
    x_edge=centre[0]+1
    dis_straight.append([x for x in range(x_edge)])
    dos_straight.append([ldos[x,y] for x in range(x_edge)])
    dos_diagonal.append([ldos[x,x] for x in range(x_edge)])
    
    dis_imp.append([(2/n_x)*x for x in range(x_edge)])
    sqrt2=np.sqrt(2)
    dis_diagonal.append([sqrt2*x for x in range(x_edge)])

########################################################
########################## Plot ########################
########################################################
def main():
    label_sites = r'No. sites from $r_0$ to $r_\text{horiz}$'
    label_diagonal = r'No. sites from $r_0$ to $r_\text{diag}$'
    label_impurity_straight = r'Distance $r_0$ to $r_\text{horiz}$'
    label_impurity_diagonal = r'Distance $r_0$ to $r_\text{diag}$'
    
    colors = [cm.rainbow(i) for i in np.linspace(0, 1, n_data)]

    fig = plt.figure()

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222, sharey=ax1)
    ax3 = fig.add_subplot(223, sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(224, sharex=ax2, sharey=ax1)

    index=sorted(range(len(n_list)), key=lambda k: n_list[k])
    
    for j,i in enumerate(index):
        n = n_list[i]
        ax1.plot(dis_straight[i], dos_straight[i], label=f'${n}$', color=colors[j])
        
        ax2.plot(dis_imp[i], dos_straight[i], label=f'${n}$', color=colors[j])

        ax3.plot(dis_diagonal[i], dos_diagonal[i], label=f'${n}$', color=colors[j])

        ax4.plot(dis_imp[i], dos_diagonal[i], label=f'${n}$', color=colors[j])
        
    handles, labels = ax1.get_legend_handles_labels()
    handles, labels = handles, labels

    lower_limit = 0 
    upper_limit = np.max([np.max(dos_straight), np.max(dos_diagonal)])*3.1
    # plt.ylim(lower_limit, upper_limit)

    ax1.set(xlabel=label_sites)
    ax2.set(xlabel=label_impurity_straight)
    ax3.set(xlabel=label_diagonal)
    ax4.set(xlabel=label_impurity_diagonal)
    # ax1.set(ylabel='LDOS')
    # ax3.set(ylabel='LDOS')
    fig.text(-0.01, 0.5, 'Local density of states', va='center', rotation='vertical')
    
    ax2.xaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax4.xaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)
    ncol = int(np.ceil(len(handles)/2))
    legend = fig.legend(handles, labels, title='Number of sites along each axis, $n\equiv n_x=n_y$', loc="upper center", mode = "expand", ncol = ncol, 
    fancybox=True, shadow=True, #prop=fontP,
    bbox_to_anchor=(0,0.95,1,0.2))

    # k_F = Fermi_vector(mu, t, omega)
    # friedel_wavelength = Friedel_wavelength(k_F)
    # text_fermi = (f'Friedel_wavelength $={friedel_wavelength:.2f}$')
    # fig.text(0.81, .14, text_fermi,
             # {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           # ec="none", pad=0.2)}, ha='center', va='center')


    # plt.subplots_adjust(wspace=-.5, hspace=-1.)

    fig.set_size_inches(w=latex_width, h=5) 

    plt.tight_layout()

    return fig
########################################################
######################### Caption ######################
########################################################
def caption():
    fermi_vector = Fermi_vector(mu, t, omega)
    friedel_wavelength = Friedel_wavelength(fermi_vector)
    [nx,ny,nz]=dimensions=data.dimensions

    text=rf'''The local density of states of a spinless square lattice tight-binding
model at $\mu/t={mu:.2f}$ with ${nx}\times{ny}$
sites, a single orbital with an impurity at the centre with coupling
strength $V/t={V:.2f}$. The impurity gives rise to Fridel's eponymous
waves in the electron quasiparticle density. The electronic excitations
at zero temperature necessarily carry the Fermi energy, and hence the 
wavefunction describing the excitations is of the Fermi wavelength. 
The electronic charge distrbution is the square modulus of the
wavefunction and hence takes on twice the periodicty or double the
wavelength $\lambda_\text{{Friedel}}=\lambda_\text{{Fermi}}/2=
{friedel_wavelength:.3}$. \\
'''
    with open(output+'.txt', 'w') as f:
        f.write(text)

##########################################################
########################## Main ##########################
########################################################## 
main()
plt.savefig(output+'.pdf', bbox_inches = "tight")
caption()
