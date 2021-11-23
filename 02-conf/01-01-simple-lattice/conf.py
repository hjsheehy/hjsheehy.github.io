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
def CreateName(n,n_spins,n_orbs):
    return f'Normal_{str(n).zfill(2)}_{n_spins:d}_{n_orbs:d}'.replace('0.00e+00','0')


def Main():
    ClearOutFolder(out_folder, silent)
    for n_orbs in n_orbs_list:
        confname = os.path.join(out_folder,CreateName(n,n_spins,n_orbs))
        temp_conf = f'''from lib import *
####################################################################
############################# Square ###############################
####################################################################
model='tb'
n_x=n_y=n={n}
n_z=1
n_spins={n_spins}
n_orbs={n_orbs}
mu={mu}
t={t}
V=1.21
'''
        with open(confname+'.conf', 'w') as f:
            f.write(temp_conf)
###################################################################
############################# Main ################################
###################################################################
n = 3
n_spins = 1
n_orbs_list = [1,2]
a = 1.0
mu= -3.57
t = 1.0
V = 0.1
impurity_locations = [[0,0]]
epsilon = 0.1
U=2.8

Main()
