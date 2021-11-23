from plt_lib import *
########################################################
########################## Data ########################
########################################################
layer=(slice(None),slice(None),0)
omega=0
orbital=0
n_data=len(filenames)

nrow = 2; ncol = 2;

fig, axs = plt.subplots(nrows=nrow, ncols=ncol)
########################################################
######################### Plot a #######################
########################################################
label_x = r'$x/a_x$'
label_y = r'$y/a_y$'

label_sites = r'No. sites from $r_0$ to $r_\text{horiz}$'
labels=['Normal state', 'BCS state', r'Equal spin n.n.', r'Spin singlet n.n.']

colors = [cm.rainbow(i) for i in np.linspace(0, 1, n_data)]
min_list,max_list=[],[]
for k in range(n_data):
    filename = filenames[k]
    globals().update(conf_file(filename))

    import _pickle as cPickle
    with open(filename, 'rb') as f:
        data = cPickle.load(f)
    
    if model=='tb':
        i=j=0
        plotdata=plot_data(fig, axs[i,j], data)
        fig, axs[i,j] = plotdata.differential_current_map(omega)
        plotdata.set_label(labels[0])
        plotdata.set_cbar()
        plotdata.set_text_box()

    if model=='sc':
        if hubbard=='onsite':
            i=0; j=1
            plotdata=plot_data(fig, axs[i,j], data)
            fig, axs[i,j] = plotdata.differential_current_map(omega)
            plotdata.set_label(labels[1])
            plotdata.set_cbar()

        if hubbard=='nn_equal':
            i=1; j=0
            plotdata=plot_data(fig, axs[i,j], data)
            fig, axs[i,j] = plotdata.differential_current_map(omega)
            plotdata.set_label(labels[2])
            plotdata.set_cbar()

        if hubbard=='nn_opposite':
            i=1; j=1
            plotdata=plot_data(fig, axs[i,j], data)
            fig, axs[i,j] = plotdata.differential_current_map(omega)
            plotdata.set_label(labels[3])
            plotdata.set_cbar()

axs[0,0].get_shared_x_axes().join(axs[0,0],axs[0,1])
axs[1,0].get_shared_x_axes().join(axs[0,0],axs[1,0])
axs[1,1].get_shared_x_axes().join(axs[1,1],axs[1,0])
axs[0,0].set_xticklabels([])
axs[0,1].set_xticklabels([])
axs[0,1].get_shared_y_axes().join(axs[0,1],axs[0,0])
axs[1,0].get_shared_y_axes().join(axs[1,1],axs[1,0])
axs[0,1].set_yticklabels([])
axs[1,1].set_yticklabels([])
axs[0,0].set(xlabel='')
axs[0,1].set(xlabel='',ylabel='')
axs[1,1].set(ylabel='')
axs[1,1].get_shared_x_axes().join(axs[0,1],axs[1,1])

plt.tight_layout()

fig.set_size_inches(w=latex_width, h=4.5) 

plt.savefig(output+'_a.pdf', bbox_inches = "tight")
plt.close()
########################################################
######################### Plot b #######################
########################################################
nrow=2; ncol=1
fig, axs = plt.subplots(nrows=nrow, ncols=ncol)

for k in range(n_data):
    filename = filenames[k]
    globals().update(conf_file(filename))

    import _pickle as cPickle
    with open(filename, 'rb') as f:
        data = cPickle.load(f)

    [nx,ny,nz]=dimensions=data.dimensions

    plotdata = plot_data(fig, axs[0], data)
    ldos = data.local_density_of_states(omega,orbital=0)
    vmin,vmax=np.amin(ldos),np.amax(ldos)
    min_list.append(vmin)
    max_list.append(vmax)

    centre=data.centre
    y=0
    x_edge=centre[0]+1
    dis_straight=[x for x in range(x_edge)]
    dos_straight=[ldos[x,y] for x in range(x_edge)]

    axs[0].plot(dis_straight, dos_straight, color=colors[k])

    if model=='tb':
        label=labels[0]
        xxx,yyy=10,-30
    if model=='sc':

        fock=data.fock()
        gorkov=data.gorkov()
        if hubbard=='onsite':
            label=labels[1]
            gorkov=np.real([gorkov[x,0,0,0,0,x,0,0,1,0] for x in range(nx)])
            fock=np.real([fock[x,0,0,0,0,x,0,0,1,0] for x in range(nx)])
            xxx,yyy=-70,-30
        if hubbard=='nn_equal':
            label=labels[2]
            gorkov=np.real([gorkov[x,0,0,0,0,(x+1)%nx,0,0,0,0] for x in range(nx)])
            fock=np.real([fock[x,0,0,0,0,(x+1)%nx,0,0,0,0] for x in range(nx)])
            xxx,yyy=-60,30
        if hubbard=='nn_opposite':
            label=labels[3]
            xxx,yyy=20,20
            gorkov=np.real([gorkov[x,0,0,0,0,(x+1)%nx,0,0,1,0] for x in range(nx)])
            fock=np.real([fock[x,0,0,0,0,(x+1)%nx,0,0,1,0] for x in range(nx)])
    
        gorkov=gorkov[:centre[0]+1]
        fock=fock[:centre[0]+1]
        axs[1].plot(dis_straight, gorkov, color=colors[k])
        axs[1].plot(dis_straight, fock, color=colors[k], linestyle='dashed')

    xx=x_edge-int(data.centre[0]/2) + xxx/40
    yy=ldos[y,x_edge]

    axs[0].annotate(label,
                xy=(xx,yy), xycoords='data',
                xytext=(xxx, yyy), textcoords='offset pixels',
        arrowprops=dict(arrowstyle="->",
                        connectionstyle="angle3,angleA=0,angleB=-90"))

# a.set_label(r'Gorkov channel')
# b.set_label(r'Fock channel')
# (axs[1]).legend()

axs[1].legend(['Gorkov channel', 'Fock channel'],labelcolor='k')

vmin,vmax=min(min_list),max(max_list)

axs[0].set(ylabel=f'LDOS($\omega=0$)'+r'$|_{{\mathbf{{r}}=[x,0]}}$')
axs[1].set(ylabel=f'Renormalisation field')
axs[1].set(xlabel=label_x)
axs[0].set_ylim([vmin,1.01*vmax])
############################################# 
############# Fock and Gorkov fields ########
#############################################
plt.tight_layout()

fig.set_size_inches(w=latex_width, h=4.5) 

plt.savefig(output+'_b.pdf', bbox_inches = "tight")
########################################################
######################### Caption ######################
########################################################
def caption():
    fermi_vector = Fermi_vector(mu, t, omega)
    friedel_wavelength = Friedel_wavelength(fermi_vector)

    text=rf'''The local density of states of a spinless 
square lattice tight-binding model at $\mu/t={mu:.2f}$ 
with ${nx}\times{ny}$ sites, a single orbital with an 
impurity at the centre with coupling strength $V/t={V:.2f}$. 
The superconducting models are solved self-consistently. 
Observe that in the BCS model (with interaction coupling 
$U={U}$), the  Friedel oscillations remain, but are 
attenuated in magnitude, as electronic charge is lost
to the formation of Cooper pairs.
In the simplest nearest-neighbour model of superconductivity,
an onsite BCS model with an additional equal-spin,
nearest-neighbour pairing term (with coupling 
$U_\text{{nn}}={U_nn}$), no equal-spin,
nearest-neighbour superconductivity emerges. Charge density 
is distributed into the Fock channel, which is a renormalisation
of the hopping parameter. The Fock channel oscillates in space
in phase with the Fridel oscillations of the charge density.
Adding a nearest-neighbour, singlet term 
($U_\text{{nn, singlet}}={U_nn_opp}$) to the previous
model, we observe nearest-neighbour Cooper pairing
and zero Fock density, like the conventional BCS model.
'''
    with open(output+'.txt', 'w') as f:
        f.write(text)
########################################################
########################## Main ########################
########################################################
caption()
