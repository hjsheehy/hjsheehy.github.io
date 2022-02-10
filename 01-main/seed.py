from lib import *
in_folder = os.path.join(MAIN,IN)

seeds=['INT','unitary']
xlabel='U_T'
ylabel='U_R'
xmin,xmax,dx=0.012,4.423,0.2002
ymin,ymax,dy=0.008,4.312,0.1982

seed_coords=seed_coord(xmin=xmin, xmax=xmax, dx=dx, ymin=ymin, ymax=ymax, dy=dy, separation=3)

create_seeds(seed_coords=seed_coords, seeds=seeds, xlabel=xlabel, ylabel=ylabel, output_folder=in_folder)
