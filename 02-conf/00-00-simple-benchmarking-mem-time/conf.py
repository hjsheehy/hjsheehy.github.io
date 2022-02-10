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
def CreateName(n,n_spins,n_orbitals):
    return f'Normal_{str(n).zfill(2)}_{n_spins:d}_{n_orbitals:d}'.replace('0.00e+00','0')

def orbs(n_orbitals):
    if n_orbitals==1:
        return '''orb0=np.array([0,0])'''
    if n_orbitals==2:
        return '''orb0=np.array([0,0])
orb1=np.array([0.35,0])'''
        return
def orbitals(n_orbitals):
    if n_orbitals==1:
        return '''orbitals=[orb0]'''
    if n_orbitals==2:
        return '''orbitals=[orb0,orb1]'''
def hop_links(n_orbitals):
    if n_orbitals==1:
        return '''hop_links = [[t, 0, 0, [1,0]],
[t, 0, 0, [0,1]]]'''
    if n_orbitals==2:
        return '''hop_links = [[t, 0, 1, [0,0]],
[t, 1, 0, [1,0]],
[t, 1, 1, [0,1]],
[t, 0, 0, [0,1]]]'''

def Main():
    ClearOutFolder(out_folder, silent)
    for n in n_list:
        for n_orbitals in n_orbitals_list:
            for n_spins in range(1,n_orbitals+1):
                confname = os.path.join(out_folder,CreateName(n,n_spins,n_orbitals))
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
{orbs(n_orbitals)}
dimensions=[n_x,n_y]
basis=[a1,a2]
######### lattice & sites #########
mu={mu}
{orbitals(n_orbitals)}
n_orbitals={n_orbitals}
SO_spin=np.eye(n_spins)
SO_orbital=np.eye(n_orbitals)
SO_tensor=-mu*np.kron(SO_spin,SO_orbital)
############ hopping ##############
t={t}
{hop_links(n_orbitals)}
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
n_list = np.arange(1,81,2)
n_spins_list = [1,2]
n_orbitals_list = [1,2]
a = 1.0
mu= -3.57
t = 1.0
V = 0.000001
impurity_locations = [[0,0]]
epsilon = 0.1

Main()
