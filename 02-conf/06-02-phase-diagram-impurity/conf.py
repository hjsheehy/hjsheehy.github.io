import sys
import os

ROOT_DIR = os.path.dirname(
            os.path.dirname(
             os.path.dirname(os.path.abspath(__file__))))
os.chdir(ROOT_DIR)
FILENAME = os.path.dirname(os.path.abspath(__file__)).split('/')[-1]
sys.path.append('01-main/')
from phase_diagram import *

INITIAL_FIELDS=['INT','unitary']
xlabel='U_v'
ylabel='U_w'
xmin,xmax,dx=0.000,4.423,0.1002
ymin,ymax,dy=0.000,4.312,0.0982
n_layers=16
separation=24
V=20.28
parent_sim='06-01*'

create_config_files(xlabel=xlabel, ylabel=ylabel, xmin=xmin, xmax=xmax, dx=dx, ymin=ymin, ymax=ymax, dy=dy, n_layers=n_layers, separation=separation, FILENAME=FILENAME, INITIAL_FIELDS=INITIAL_FIELDS, V=V, parent_sim=parent_sim)
