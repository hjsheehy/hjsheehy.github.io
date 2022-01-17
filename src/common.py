'''Common functions for large projects'''
import os
import string
import numpy as np
import random

def GetCh():
    """Read single character from standard input without echo."""
    import sys, tty, termios
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        ch = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    return ch

def YesNo(question):
    c = ""
    print(question + " [Y/n]: ", end = "", flush = True)
    while c not in ("y", "Y", "n", "N"):
        c = GetCh().lower()
    return c == 'y' or c == "Y"

def CreateFolder(folder):
    '''If folder does not exist, create one.'''
    if os.path.exists(folder):
        pass
    else:
        print(f'{folder}\nNo folder found. Created one.')
        os.system(f'mkdir {folder}')

def ClearOutFolder(out_folder, silent=False):
    '''If out_folder does not exists, creates one. Prompts clearing
    of out_folder contents if not-empty.'''
    CreateFolder(out_folder)
    if len(os.listdir(out_folder)) != 0:
        if not silent:
            if YesNo(f'Clear out folder before creating new .conf files located at\n{out_folder}?'):
                pass
            else:
                exit()
        os.system(f'rm -r {out_folder}*')

suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']

def humansize(nbytes):
    '''Returns human readable memory'''
    i = 0
    while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024.
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])

def closest_point(points,point):
    points = np.asarray(points)
    deltas = points - point
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)

def create_seed_coords(xmin, xmax, dx, ymin, ymax, dy, separation=7):
    """Creates n_seeds random 2D coordinates on the grid defined by the bounds, with a minimum separation. e.g. to seed inital points for Monte Carlo."""
    coords=[]
    x=np.arange(xmin,xmax,dx)
    y=np.arange(ymin,ymax,dy)
    x=x[::separation]
    y=y[::separation]
    for xx in x:
        for yy in y:
            coords.append([xx,yy])
    return coords

def create_layer_coords(seed_coord, layer_index, xmin, xmax, dx, ymin, ymax, dy, separation=7):
    """Creates n_seeds random 2D coordinates on the grid defined by the bounds, with a minimum separation. e.g. to seed inital points for Monte Carlo."""
    i=layer_index
    coords=np.indices([2*i+1,2*i+1])-i
    coords[:,1:-1,1:-1]=0
    coords=coords.T
    coords=coords.reshape(-1, coords.shape[-1])
    coords = np.array([coord for coord in coords if not np.array_equal(coord,np.array([0,0]))])
    [imin,imax]=[-separation,separation]
    if len(coords)==0:
        return []
    coords = np.array([coord for coord in coords if (imin<coord[0]<imax and imin<coord[1]<imax)])
    if len(coords)==0:
        return []
    coords=np.multiply(coords,[dx,dy])
    coords=coords+np.array(seed_coord)
    coords = np.array([coord for coord in coords if (xmin-dx<coord[0]<xmax+dx and ymin-dy<coord[1]<ymax+dy)])
    return coords

def create_seed(seed_index, seed_coord, initial_fields_label, xlabel, ylabel, output_folder):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    temp_conf = f'''initial_fields_label='{initial_fields_label}'
[{xlabel}, {ylabel}]={seed_coord}
'''
    confname=f'{0:03d}_{seed_index:03d}_'
    for coord in seed_coord:
        confname+=f'{coord:.3f}_'
    confname+=initial_fields_label
    confname = os.path.join(output_folder,confname)
    with open(confname+'.conf', 'w') as f:
        f.write(temp_conf)

def create_seeds(seed_coords, initial_fields_labels, xlabel, ylabel, output_folder):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for initial_fields_label in initial_fields_labels:
        for seed_index,seed_coord in enumerate(seed_coords):

            create_seed(seed_index=seed_index, seed_coord=seed_coord, initial_fields_label=initial_fields_label, xlabel=xlabel, ylabel=ylabel, output_folder=output_folder)

def create_germ(layer_coord, dx, dy, seed_index, seed_coord, layer_index, xlabel, ylabel, output_folder):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    parent_filename=parent_label(seed_coord, seed_index, layer_coord, layer_index, dx, dy)
    temp_conf = f'''parent_filename='{parent_filename}'
[{xlabel}, {ylabel}]={list(layer_coord)}
'''
    confname=f'{layer_index:03d}_{seed_index:03d}'
    for coord in layer_coord:
        confname+=f'_{coord:.3f}'
    confname = os.path.join(output_folder,confname)
    with open(confname+'.conf', 'w') as f:
        f.write(temp_conf)

def create_germs(germ_coords, dx, dy, seed_index, seed_coord, layer_index, xlabel, ylabel, output_folder):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for layer_coord in germ_coords:
        create_germ(layer_coord=layer_coord, dx=dx, dy=dy, seed_index=seed_index, seed_coord=seed_coord, layer_index=layer_index, xlabel=xlabel, ylabel=ylabel, output_folder=output_folder)

def germ_parent_coord(seed_coord, germ_coord, dx, dy):
    a=np.around(np.sqrt(1/2),4)
    n=np.array(germ_coord-seed_coord)/np.array([dx,dy])
    n=np.around(n/np.linalg.norm(n),4)
    Px=np.dot(n,[1,0])
    Py=np.dot(n,[0,1])
    Px=np.heaviside(np.abs(Px)-a,1)*dx*np.sign(Px)
    Py=np.heaviside(np.abs(Py)-a,1)*dy*np.sign(Py)
    parents_coord=np.array(germ_coord)-np.array([Px,Py])
    parents_coord=np.around(parents_coord,3)
    return parents_coord

xmin,xmax,dx=0.012,4.423,0.2002
ymin,ymax,dy=0.008,4.312,0.1982
seed_coord=np.array([0.012,1.395])
germ_coord=np.array([0.212,1.197])
parent_coord=germ_parent_coord(seed_coord, germ_coord, dx, dy)

def parent_label(seed_coord, seed_index, germ_coord, germ_layer_index, dx, dy):
    parent_layer=germ_layer_index-1
    confname=f'{parent_layer:03d}_{seed_index:03d}_'
    parent_coord=germ_parent_coord(seed_coord=seed_coord, germ_coord=germ_coord, dx=dx, dy=dy)
    for coord in parent_coord:
        confname+=f'{coord:.3f}_'
    confname=confname+'*'
    return confname

def layer_no(confName):
    bud=confName.split('/')[-1]
    bud=bud.split('_')[0]
    return int(bud)

def is_seed(confName):
    return layer_no(confName)==0

def active_label(label):
    return 'active_'+label

def inactive_label(label):
    return 'inactive_'+label.split('active_')[1]
