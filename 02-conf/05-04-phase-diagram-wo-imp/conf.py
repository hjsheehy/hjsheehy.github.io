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
####################################################################
########################## Configuration ############################
####################################################################
def CreateName(mu,U_R,U_T,n):
    return f'lattice_{mu:0.2f}_{U_R:0.4f}_{U_T:0.4f}_{n:d}'

def Main():
    ClearOutFolder(out_folder)
    for mu in MU:
        for U_T in T:
            for U_R in R:
                confname = os.path.join(out_folder,CreateName(mu,U_R,U_T,n))
                temp_conf = f'''n={n}
mu={mu}
V=2.21
w=0.79
v=0
U_T={U_T}
U_R={U_R}
phi=4.32
zeta_v=4.52
zeta_w=3.23
Delta_v=5.92
Delta_w=4.73'''
                with open(confname+'.conf', 'w') as f:
                    f.write(temp_conf)
###################################################################
############################# Main ################################
###################################################################
n=41
#MU=[-3.92,-2.43,-1.32,-0.32,0.00,0.32,1.32,2.43,3.92]
MU=[-3.92,-0.32,0.32,0.00,2.12,3.92]
R=np.arange(0.012,4.423,0.1002)
T=np.arange(0.008,4.312,0.0982)
Main()
