from plt_lib import *
########################################################
########################## Data ########################
########################################################
layer=(slice(None),slice(None),0)
omega=0
n_data=len(filenames)

nrow = 2; ncol = 2;

fig, axs = plt.subplots(nrows=nrow, ncols=ncol)

# ldos=[]
# for i in range(n_data):
#     filename = filenames[i]
#     globals().update(conf_file(filename))

#     import _pickle as cPickle
#     with open(filename, 'rb') as f:
#         model = cPickle.load(f)
########################################################
########################## Plot ########################
########################################################
def main():
    label_x = r'$x/a_x$'
    label_y = r'$y/a_y$'
    
    label_sites = r'No. sites from $r_0$ to $r_\text{horiz}$'
    labels=['TB', 'BCS', 'BCS n.n.']

    colors = [cm.rainbow(i) for i in np.linspace(0, 1, n_data)]
    min_list,max_list=[],[]

    for k in range(n_data):
        filename = filenames[k]
        globals().update(conf_file(filename))


        import _pickle as cPickle
        with open(filename, 'rb') as f:
            Model = cPickle.load(f)
        ldos = LDOS(Model.density_of_states, Model.omegas, omega, trace_over=True, layer=layer)
        vmin,vmax=np.amin(ldos),np.amax(ldos)
        min_list.append(vmin)
        max_list.append(vmax)

        centre=Model.centre
        y=0
        x_edge=centre[0]+1
        dis_straight=[x for x in range(x_edge)]
        dos_straight=[ldos[x,y] for x in range(x_edge)]

        i=j=1
        axs[i,j].plot(dis_straight, dos_straight, label=f'$n$', color=colors[k])
        
        xx=x_edge-5
        yy=ldos[y,x_edge]

        if model=='tb':
            label=labels[0]
            xxx,yyy=-40,-20
        if model=='sc':
            if hubbard=='onsite':
                label=labels[1]
                xxx,yyy=-40,-20
            if hubbard=='nn':
                label=labels[2]
                xxx,yyy=-50,20

        axs[i,j].annotate(label,
                    xy=(xx,yy), xycoords='data',
                    xytext=(xxx, yyy), textcoords='offset pixels',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="angle3,angleA=0,angleB=-90"))

        if model=='tb':
            i=j=0
            ldos_tb = LocalDensityOfStates(fig, axs[i,j], Model, ldos=ldos, omega=omega)
        if model=='sc':
            if hubbard=='onsite':
                i=0; j=1
                ldos_bcs = LocalDensityOfStates(fig, axs[i,j], Model, ldos=ldos, omega=omega)
            if hubbard=='nn':
                i=1; j=0
                ldos_nn = LocalDensityOfStates(fig, axs[i,j], Model, ldos=ldos, omega=omega)

    vmin,vmax=min(min_list),max(max_list)
    
    for i,x in enumerate([ldos_tb,ldos_bcs,ldos_nn]):
        # x.vmin=vmin
        # x.vmax=vmax
        x.imshow()
        x.set_cbar()
        x.set_label(labels[i])
    
    ldos_tb.set_text_box()

    axs[0,0].get_shared_x_axes().join(axs[0,0],axs[0,1])
    axs[0,0].set_xticklabels([])
    axs[0,1].get_shared_y_axes().join(axs[0,1],axs[0,0])
    axs[0,1].set_yticklabels([])
    axs[0,0].set(xlabel='')
    axs[0,1].set(xlabel='',ylabel='')
    axs[1,1].set(xlabel=label_x, ylabel=f'\n\n\nLDOS($\omega=0$)'+r'$|_{{\mathbf{{r}}=[x,0]}}$')
    #axs[1,1].tick_params(axis="y",direction="in", pad=-22)
    #axs[1,1].get_shared_x_axes().join(axs[0,1],axs[1,1])
    axs[1,1].set_ylim([0,1.01*vmax])

    plt.tight_layout()

    fig.set_size_inches(w=latex_width, h=4.5) 
    return
########################################################
######################### Caption ######################
########################################################
def caption():
    fermi_vector = Fermi_vector(mu, t, omega)
    friedel_wavelength = Friedel_wavelength(fermi_vector)

    text=rf'''The local density of states of a spinless square lattice tight-binding
model at $\mu/t={mu:.2f}$ with ${dimensions[0]}\times{dimensions[1]}$
sites, a single orbital with an impurity at the centre with coupling
strength $V/t={V:.2f}$. The impurity gives rise to Fridel's eponymous
waves in the electron quasiparticle density. The electronic excitations
at zero temperature necessarily carry the Fermi energy, and hence the 
wavefunction describing the excitations is of the Fermi wavelength. 
The electronic charge distrbution is the square modulus of the
wavefunction and hence takes on twice the periodicty or double the
wavelength $\lambda_\text{{Friedel}}=\lambda_\text{{Fermi}}/2=
{friedel_wavelength:.3}$.
'''
    with open(output+'.txt', 'w') as f:
        f.write(text)
########################################################
########################## Main ########################
########################################################
main()
plt.savefig(output+'.pdf', bbox_inches = "tight")
caption()
