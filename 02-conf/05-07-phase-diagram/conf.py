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
silent=sys.argv[-1]
ClearOutFolder(out_folder,silent=silent)

initial_fields_labels=['INT','unitary']
xlabel='U_T'
ylabel='U_R'
xmin,xmax,dx=0.012,4.423,0.2002
ymin,ymax,dy=0.008,4.312,0.1982
n_layers=5
separation=7

seed_coords=create_seed_coords(xmin=xmin, xmax=xmax, dx=dx, ymin=ymin, ymax=ymax, dy=dy, separation=separation)

create_seeds(seed_coords=seed_coords, initial_fields_labels=initial_fields_labels, xlabel=xlabel, ylabel=ylabel, output_folder=out_folder)

tmp=int(separation/2)
if n_layers<tmp:
    raise ValueError("n_layers must be greater than half the separation distance! Using n_layers=separation/2 instead.")
    n_layers=tmp

seed_coord=[1,0]
layer_index=1
layer_coords=create_layer_coords(seed_coord=seed_coord, layer_index=layer_index, xmin=xmin, xmax=xmax, dx=dx, ymin=ymin, ymax=ymax, dy=dy, separation=separation)

for seed_index,seed_coord in enumerate(seed_coords):
    for layer_index in range(1,n_layers):
        germ_coords=create_layer_coords(seed_coord=seed_coord, layer_index=layer_index, xmin=xmin, xmax=xmax, dx=dx, ymin=ymin, ymax=ymax, dy=dy, separation=separation)
        create_germs(germ_coords=germ_coords, dx=dx, dy=dy, seed_index=seed_index, seed_coord=seed_coord, layer_index=layer_index, xlabel=xlabel, ylabel=ylabel, output_folder=out_folder)
