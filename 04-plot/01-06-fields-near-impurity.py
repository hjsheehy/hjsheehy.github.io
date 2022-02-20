from plt_lib import *
########################################################
########################## Data ########################
########################################################
layer=(slice(None),slice(None),0)
omega=0
n_data=len(filenames)

n_x_tb=[]
n_x_sc=[]
n_z_tb=[]
n_z_sc=[]
ldos_x_tb=[]
ldos_z_tb=[]
gorkov_x_sc=[]
gorkov_z_sc=[]

for i in range(n_data):
    ################ import ################
    filename = filenames[i]
    globals().update(conf_file(filename))

    import _pickle as cPickle
    with open(filename, 'rb') as f:
        Model = cPickle.load(f)
    ldos = LDOS(Model.density_of_states, Model.omegas, omega, trace_over=True, layer=layer)
    
    y=z=0
    m=0
    indices_up=Model.index[:int(n_x/2),y,z,0,m]
    indices_down=Model.index[:int(n_x/2),y,z,1,m]
    print(model)
    print(Model.hartree[indices_up])
    print(Model.hartree[indices_down])
    print(np.shape(Model.fock))
    exit()

    if model=='tb':
        if n_name=='n':
            n_x_tb.append(n_x)
            ldos_x_tb.append(ldos[0,0])
        if n_name=='n_z':
            n_z_tb.append(n_z)
            ldos_z_tb.append(ldos[0,0])
    if model=='sc':
        if n_name=='n':
            n_x_sc.append(n_x)
            gorkov_x_sc.append(Model.gorkov[0])
        if n_name=='n_z':
            n_z_sc.append(n_z)
            gorkov_z_sc.append(Model.gorkov[0])
            n_fixed=n_x

nrow = 2; ncol = 1;

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
    label_x = r'$n_x=n_y$, with $n_z=1$'
    label_z = f'$n_z$, with $n_x=n_y={n_fixed}$'
    label_y0 = f'LDOS($\omega=0$)\nTightbinding model'
    label_y1 = r'$\Delta_{{\uparrow,\downarrow}}/t$'+f'\nSelf-consistent BCS'
    label_y2= label_y0
    label_y3= label_y1

    x0=np.array(n_x_tb)
    y0=np.array(ldos_x_tb)
    x1=np.array(n_x_sc)
    y1=np.array(gorkov_x_sc)
    x2=np.array(n_z_tb)
    y2=np.array(ldos_z_tb)
    x3=np.array(n_z_sc)
    y3=np.array(gorkov_z_sc)
    
    inds = x0.argsort()
    x0=x0[inds]
    y0=y0[inds]
    inds = x1.argsort()
    x1=x1[inds]
    y1=y1[inds]
    inds = x2.argsort()
    x2=x2[inds]
    y2=y2[inds]
    inds = x3.argsort()
    x3=x3[inds]
    y3=y3[inds]

    print(y2)
    print(y3)

    i=0
    
    fig = plt.figure()

    color='tab:red'
    ax0 = fig.add_subplot(211)
    ax0.plot(x0, y0,color=color,marker='.',markersize=4,label=label_y0)
    ax0.set_ylabel(label_y0, color=color)
    ax0.tick_params(axis='y', labelcolor=color)
    ax0.set_xlabel(label_x)

    ax1 = ax0.twinx()

    color = 'tab:blue'
    ax1.plot(x1, y1,color=color,marker='x',markersize=4,label=label_y1)  
    ax1.set_ylabel(label_y1, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    color='tab:red'
    ax2 = fig.add_subplot(212)
    ax2.plot(x2, y2,color=color,marker='.',markersize=4,label=label_y2)
    ax2.set_ylabel(label_y2,color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_xlabel(label_z)
    
    ax3 = ax2.twinx()

    color = 'tab:blue'
    ax3.plot(x3, y3,color=color,marker='x',markersize=4,label=label_y3)  
    ax3.tick_params(axis='y', labelcolor=color)
    ax3.set_ylabel(label_y3, color=color)

    ax0.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()

    fig.set_size_inches(w=LATEX_WIDTH, h=4.5) 
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
