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
def CreateName(*int_list):
    name=""
    for a in list(int_list):
        name = name + f'{a:03d}_'
    return name[:-1]

def Main():
    ClearOutFolder(out_folder)
    N=len(T)
    for i0,state in enumerate(states):
        for i1,V in enumerate(VV):
            for i2,mu in enumerate(MU):
                for i3,trajectory in enumerate(trajectories):
                    for i4,U_R in enumerate(R):
                        for i5,U_T in enumerate(T):
                            if trajectory=='forward':
                                pass
                            elif trajectory=='reverse':
                                i5=N-i5-1
                            confname = os.path.join(out_folder,CreateName(i0,i1,i2,i3,i4,i5))
                            temp_conf = f'''n={n}
mu={mu}
V={V}
w=0.79
v=1.31
U_T={U_T}
U_R={U_R}
phi=6.32
zeta_v=1.52
zeta_w=1.23
Delta_v=3.52
Delta_w=2.23
state='{state}'
trajectory='{trajectory}'
'''
                            with open(confname+'.conf', 'w') as f:
                                f.write(temp_conf)
###################################################################
############################# Main ################################
###################################################################
n=41
#MU=[-3.92,-2.43,-1.32,-0.32,0.00,0.32,1.32,2.43,3.92]
# MU=[-3.92,-0.32,0.32,0.00,2.12,3.92]
MU=[-2.12]
R=[4.423,3.021,1.121]
T=[4.312,1.471,3.311]
VV=[1.47]
states=['INT']
trajectories=['reverse']

Main()
