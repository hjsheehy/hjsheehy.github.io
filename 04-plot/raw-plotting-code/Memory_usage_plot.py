import os
import numpy as np 
import matplotlib
from matplotlib import pyplot as plt
import scipy as sp
from scipy.special import j0, y0
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})

name = os.path.basename(__file__)
name = (name.split('.'))[0]

def Main():
    x = np.linspace(0, 121, 200) 
    def f(N):
        return 16*(2**3*N**2)**2*10**(-9)
    # 16 bits per complex data type 
    # 2 spin, 2 orbital, 2 BdG, squared since matrix
    # 10^(-9) Gigabytes

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, f(x)) 

    fig.set_size_inches(w=latex_width, h=3.5) 
    ax.set(xlabel=r"Number of sites, N", ylabel=r"Memory, gb")
    plt.savefig('out/'+name + '.pgf')
###############################################
#################### Main #####################
###############################################
latex_width=4.7747
Main()