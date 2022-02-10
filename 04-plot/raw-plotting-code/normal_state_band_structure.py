import os
import numpy as np 
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import MultipleLocator 
plt.rcParams['text.usetex'] = True
###############################################
#### 3D plots aren't well supported by pdg ####
############## Using pdf instead ##############
###############################################
# name.py outputs name.pdf
name = os.path.basename(__file__)
name = (name.split('.'))[0]
def Main():
    #######################################
    ################# Data ################
    #######################################
    def epsilon(x,y):
        return -2*(np.cos(np.pi*x)+np.cos(np.pi*y))
    X = Y = np.linspace(-1, 1, 101)
    X, Y = np.meshgrid(X, Y)
    Z = epsilon(X,Y)
    def FS(x,y,mu):
        return np.heaviside(-2*(np.cos(np.pi*x)+np.cos(np.pi*y))-mu,1)
    M0 = FS(X,Y,mu0)
    M1 = FS(X,Y,mu1)
    def QPI(x,y):
        return np.maximum(np.abs(x),np.abs(y))
    M0 = FS(X,Y,mu0)
    M1 = FS(X,Y,mu1)
    M2 = QPI(X,Y)
    #######################################
    ############# Plot data ###############
    #######################################
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    #3D:
    ax.plot_surface(X, Y, Z, cmap="jet", lw=0, rstride=1, cstride=1, linewidth=0)
    #2D projection:
    offset=-6
    ax.contourf(X, Y, Z, 50, cmap="jet", offset=offset)
    ax.contour(X, Y, Z, 10, colors="k", linewidths=0.2, linestyles="solid", offset=offset)
    ax.contour(X, Y, M0, 1, colors="white", linewidths=1, linestyles="solid", offset=offset)
    ax.contour(X, Y, M1, 1, colors="pink", linewidths=1, linestyles="solid", offset=offset)
    ax.contour(X, Y, M2, 1, colors="black", linewidths=1.5, linestyles="dashed", offset=offset)
    #######################################
    ######## Plot customisation ###########
    #######################################
    #width w from latex command \printinunitsof{in}\prntlen{\textwidth}
    fig.set_size_inches(w=latex_width, h=latex_width) 
    ax.set(#title='Band structure',
    xlabel=r"$k_x a$", ylabel=r"$k_y a$")
    # Text label bug: hidden by contour
    # ax.text(-1, 1, -7, 'FS', color='black')
    # tick labels
    label_loc=[-1,0,1]
    labels=[r'$-\pi$','0',r'$\pi$']
    plt.xticks(label_loc, labels)
    plt.yticks(label_loc, labels)
    #ax.set(zlabel=r'$\epsilon(\mathbf{k})/t$')
    ax.text(X.min()*1.1, Y.min()*1.2, Z.max()*1.3, r'$\epsilon_0(\mathbf{k})/t$')
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(offset, 4)
    # Prevent label rotation
    ax.xaxis.set_rotate_label(False)
    ax.yaxis.set_rotate_label(False)
    ax.zaxis.set_rotate_label(False)
    # Clear grid & create an empty cubic wireframe 
    ax.grid(False)
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    # remove spacing between tick labels and axes3d
    [t.set_va('center') for t in ax.get_yticklabels()]
    [t.set_ha('left') for t in ax.get_yticklabels()]
    [t.set_va('center') for t in ax.get_xticklabels()]
    [t.set_ha('right') for t in ax.get_xticklabels()]
    [t.set_va('center') for t in ax.get_zticklabels()]
    [t.set_ha('left') for t in ax.get_zticklabels()]
    # ticks pointing into the plot that leave a clean edge around the boundary
    ax.xaxis._axinfo['tick']['inward_factor'] = 0
    ax.xaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.yaxis._axinfo['tick']['inward_factor'] = 0
    ax.yaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.zaxis._axinfo['tick']['inward_factor'] = 0
    ax.zaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.zaxis._axinfo['tick']['outward_factor'] = 0.4
    #
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.zaxis.set_major_locator(MultipleLocator(4))
    ax.view_init(elev=25, azim=130)
    # fig.tight_layout()
    return fig, ax
###############################################
#################### Main #####################
###############################################
mu0=-3.57
mu1=-2.67
latex_width=4.7747
fig, ax = Main()
plt.savefig('out/'+name + '.pdf')
# plt.show()