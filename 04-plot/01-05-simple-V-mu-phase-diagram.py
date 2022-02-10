from plt_lib import *
import scipy.interpolate
########################################################
########################## Data ########################
########################################################
layer=0
omega=0
n_data=len(filenames)
y=[]
x=[]
z=[]

for i in range(n_data):
    ################ import ################
    filename = filenames[i]
    globals().update(conf_file(filename))

    data = np.load(filename, allow_pickle=True)
    dos, exec_time, mem = data['dos'], data['exec_time'], data['mem']
    ldos = LDOS(dos, omegas, omega, trace_over=True, layer=layer)
    
    y.append(V)
    x.append(mu)
    z.append(ldos[0,0])
x=np.array(x)
y=np.array(y)
z=np.array(z)
########################################################
########################## Plot ########################
########################################################
def main():
    title=r'$-\frac{1}{\pi}\Im\hat{G}^R(\omega+i\epsilon;\mathbf{r},\mathbf{r})$'+r'$|_{{\omega={}, \mathbf{{r}}=[0,0]}}$'.format(omega)
    xlabel=r'$\mu/t$'
    ylabel=r'$V/t$'
    # Set up a regular grid of interpolation points
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)

    # Interpolate
    rbf = scipy.interpolate.Rbf(x, y, z, function='gaussian')
    zi = rbf(xi, yi)
    
    fig, ax = plt.subplots(1,1,figsize=(LATEX_WIDTH,2.5))
    img = ax.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
               extent=[x.min(), x.max(), y.min(), y.max()])
    # ax.scatter(x, y, c=z, s=4, marker='s')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.title(title)
    # fig.colorbar(img)
    return fig
########################################################
######################### Caption ######################
########################################################
def caption():
    fermi_vector = Fermi_vector(mu, t, omega)
    friedel_wavelength = Friedel_wavelength(fermi_vector)

    text=rf'''
'''
    with open(output+'.txt', 'w') as f:
        f.write(text)

##########################################################
########################## Main ##########################
########################################################## 
main()
plt.savefig(output+'.pdf', bbox_inches = "tight")
caption()
