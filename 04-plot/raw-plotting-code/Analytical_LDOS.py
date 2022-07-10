from lib_plt import *

name = os.path.basename(__file__)
name = (name.split('.'))[0]

def Main():
    x = np.linspace(0, 25, 200) 
    y = [Analytical_Friedel(xx, V, mu, t, omega) for xx in x]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y,color='r') 

    fig.set_size_inches(w=latex_width, h=3.5) 
    # ax.set(title=r"LDOS$(\mathbf{r}, \omega)", xlabel=r"$r_x/a$", ylim=[0.073,0.081])
    plt.tight_layout()
    plt.savefig('out/' +name + '.pdf')
###############################################
#################### Main #####################
###############################################
latex_width=4.7747
mu=-3.88
t=1.0
omega=0.0
V=0.1
Main()
