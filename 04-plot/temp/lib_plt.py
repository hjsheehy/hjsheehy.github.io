import sys
# sys.path.append('../master/main')
from lib import *
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
# from matplotlib.font_manager import FontProperties
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib
import os
# plt.rcParams['text.usetex'] = True
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"
# matplotlib.use("pgf")
# matplotlib.rcParams.update({
    # "pgf.texsystem": "pdflatex",
    # 'font.family': 'serif',
    # 'text.usetex': True,
    # 'pgf.rcfonts': False,
# })
from scipy.special import j0, y0
from scipy.special import hankel1
from math import log10, floor, ceil
#
def round_down(x):
    return round(x, -int(floor(log10(abs(x)))))
#
def round_up(x):
    return round(x, -int(ceil(log10(abs(x)))))
#
def Analytical_memory(n, n_spins, n_orbitals):
    '''Memory (gb) of complex (16 byte) matrix representation of BdG Hamiltonian
    for system with n_spin, n_orbital and n sites along axis.'''
    return 16*np.square(np.square(n)*n_spins*n_orbitals*2)*10**(-9)
#
def Fermi_vector(mu, t, omega):
    return np.sqrt((mu+omega)/t+4)
#
def Friedel_wavelength(k_F):
    return np.pi/k_F
#   
def Analytical_Friedel(r, V, mu, t, omega):
    m_eff = 1/(2*t)
    k_omega = Fermi_vector(mu, t, np.real(omega))
    rho = m_eff/(2*np.pi)
    const = V*(m_eff**2)/(2*np.pi)
    Bessel_J = j0(k_omega*la.norm(r))
    Bessel_Y = y0(k_omega*la.norm(r))
    pert = const*Bessel_J*Bessel_Y
    return rho + pert
#
def Free_Greens_Function(r, mu, t, omega):
    m_eff = 1/(2*t)
    k_omega = Fermi_vector(mu, t, omega)
    return -1.0j*(m_eff/2)*hankel1(0, k_omega*la.norm(r))
#
def Analytical_LDOS(r, V, mu, t, omega):
    zero = np.finfo(float).eps
    rho = Free_Greens_Function(zero, mu, t, omega)
    pert = V*(Free_Greens_Function(r, mu, t, omega))**2
    return -(1/np.pi)*np.imag(rho + pert)
#
def Group_array(array,group_axis,ordering_axis):
    '''Groups array in ordering of specified group_axis'''
    index = np.argsort(array)[group_axis]
    temp = array[:,index]
    uu=np.unique(temp[group_axis,:], return_index=True)[1]
    temp = np.array([np.split(temp[i,:], uu)[1:] for i in range(np.shape(temp)[0])])
    for i in range(np.shape(temp)[ordering_axis]):
        vals = temp[ordering_axis,i,:]
        sort_index = np.unravel_index(np.argsort(vals), vals.shape)
        temp[:,i,:] = temp[:,i,sort_index][:,0,:]
    return temp
#
def ildos_line_plt_magnification(numpy_vars,confname):
    label_sites = r'Distance from impurity (sites)'
    label_impurity_straight = r'Distance to nearest neighbour'
    label_impurity_diagonal = r'Distance to next-nearest neighbour'

    data_vals = list(numpy_vars.values())
    data_keys = list(numpy_vars.keys())
    data_length = len(data_keys)

    omega=complex((data_keys[0].split('_')[-1])[0:-4]) 

    fig = plt.figure()

    ax = fig.add_subplot(111)
    for i in range(data_length): 
        data = data_vals[i]
        data_label = data_keys[i].split('_')[0]
        straight_x_sites = data['straight_x_sites']
        straight_y_DOS = data['straight_y_DOS']
        ax.plot(straight_x_sites, straight_y_DOS, label=data_label)
    
    ax.set(xlabel=label_sites)
    ax.set(ylabel='Path to nearest neighbour')
    
    handles, labels = ax.get_legend_handles_labels()
    ncol = int(np.ceil(len(handles)/2))

    inner_ax = fig.add_axes([0.59, 0.53, 0.35, 0.35])
    i=-1
    data = data_vals[i]
    data_label = data_keys[i].split('_')[0]
    straight_x_sites = data['straight_x_sites']
    straight_y_DOS = data['straight_y_DOS']
    inner_ax.plot(straight_x_sites, straight_y_DOS, label=data_label, c='y')
    inner_ax.set(title=f'$N=43$')
        
    # k_F = Fermi_vector(mu, t, omega)
    # Friedel_wavelength = Friedel_wavelength(k_F)
    # text_fermi = (f'Friedel_wavelength $={Friedel_wavelength:.2f}$')
    # fig.text(0.81, .14, text_fermi,
             # {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           # ec="none", pad=0.2)}, ha='center', va='center')

    legend = fig.legend(handles, labels, loc="upper center", mode = "expand", ncol = ncol, 
    fancybox=True, shadow=True, #prop=fontP,
    bbox_to_anchor=(0,0.93,1,0.2))


    # plt.subplots_adjust(wspace=-.5, hspace=-1.)

    fig.set_size_inches(w=latex_width, h=4) 

    plt.tight_layout()
    return fig,omega
#
latex_width=4.7747
