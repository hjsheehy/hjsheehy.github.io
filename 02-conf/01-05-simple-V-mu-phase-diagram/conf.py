import sys
import os
ROOT_DIR = os.path.dirname(
            os.path.dirname(
             os.path.dirname(os.path.abspath(__file__))))
os.chdir(ROOT_DIR)
sys.path.append('01-main/')
from lib import *
fileName = os.path.dirname(os.path.realpath(__file__)).split('/')[-1]
out_folder = os.path.join(ROOT_DIR,CONF,fileName,OUT)

silent=False
if sys.argv[-1]=='-s':
    silent=True
####################################################################
########################## Configuration ############################
####################################################################
def CreateName(n,mu,t,V):
        return f'Normal_{n:d}_{mu:0.2f}_{t:0.2f}_{V:0.2e}'.replace('0.00e+00','0')

def Main():
    ClearOutFolder(out_folder, silent)
    for mu in mu_list:
        for V in V_list:
            confname = os.path.join(out_folder,CreateName(n,mu,t,V))
            temp_conf = f'''from lib import *
####################################################################
############################# Square ###############################
####################################################################
model='tb'
n_spins={n_spins}
n_x=n_y={n}
a={a}

a1=a*np.array([1,0])
a2=a*np.array([0,1])
zero=np.array([0,0])
dimensions=[n_x,n_y]
basis=[a1,a2]
######### lattice & sites #########
mu={mu}
orbitals=[zero]
n_orbitals=len(orbitals)
SO_spin=np.eye(n_spins)
SO_orbital=np.eye(n_orbitals)
SO_tensor=-mu*np.kron(SO_spin,SO_orbital)
############ hopping ##############
t={t}
hop_links = [[t, 0, 0, [1,0]],
             [t, 0, 0, [0,1]]]
########### impurities #############
V={V}
impurity_loc = {impurity_locations}
impurity_spin = np.arange(n_spins)
impurity_orb = np.arange(n_orbitals)
impurities = [V, impurity_loc, impurity_spin, impurity_orb]
############## Energy ###############
epsilon={epsilon}
omega_lim=10*t 
omegas=np.array([0],dtype=complex)+1.0j*epsilon'''

            with open(confname+'.conf', 'w') as f:
                f.write(temp_conf)
###################################################################
############################# Main ################################
###################################################################
n_spins = 1
n = 55
a = 1.0
mu_list = np.linspace(-4.2,4.2,50)
t = 1.0
V_list = np.linspace(-2.2,2.2,50)
impurity_locations = [[0,0]]
epsilon = 0.1

Main()
