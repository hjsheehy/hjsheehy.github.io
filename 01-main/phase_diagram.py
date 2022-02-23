from lib import *
#######################################################
def save_data(bdg_model, CONFNAME, SIM_NAME):
    """Saves the data to the DATA folder, and renames to include the calculated free energy"""

    oldCONFNAME=os.path.basename(CONFNAME)+'.conf'

    free_energy=bdg_model.free_energy[-1]
    if not bdg_model.converged:
        free_energy+=100
    CONFNAME=CONFNAME+f'_{free_energy}'

    with open(CONFNAME+'.npz', 'wb') as f:
        cPickle.dump(bdg_model, f)

    CONFNAME=os.path.basename(CONFNAME)+'.conf'
    SIM_NAME=os.path.basename(SIM_NAME)
    inConf=os.path.join(MAIN,ACTIVE,SIM_NAME,oldCONFNAME)
    outConf=os.path.join(CONF,SIM_NAME,OUT,CONFNAME)
    os.rename(inConf,outConf)

def create_config_files(INITIAL_FIELDS, xlabel, ylabel, xmin, xmax, dx, ymin, ymax, dy, n_layers, separation, FILENAME, V, parent_sim):
    """Creates configuration files for the phase diagram"""

    out_folder = os.path.join(ROOT_DIR,CONF,FILENAME,OUT)
    silent=sys.argv[-1]
    ClearOutFolder(out_folder,silent=silent)

    seed_coords=create_seed_coords(xmin=xmin, xmax=xmax, dx=dx, ymin=ymin, ymax=ymax, dy=dy, separation=separation)

    create_seeds(seed_coords=seed_coords, INITIAL_FIELDS=INITIAL_FIELDS, xlabel=xlabel, ylabel=ylabel, output_folder=out_folder, V=V, parent_sim=parent_sim)

    tmp=int(separation/2)
    if n_layers<tmp:
        raise ValueError("n_layers must be greater than half the separation distance! Using n_layers=separation/2 instead.")
        n_layers=tmp

    for seed_index,seed_coord in enumerate(seed_coords):
        for layer_index in range(1,n_layers):
            germ_coords=create_layer_coords(seed_coord=seed_coord, layer_index=layer_index, xmin=xmin, xmax=xmax, dx=dx, ymin=ymin, ymax=ymax, dy=dy, separation=separation)
            create_germs(germ_coords=germ_coords, dx=dx, dy=dy, seed_index=seed_index, seed_coord=seed_coord, layer_index=layer_index, xlabel=xlabel, ylabel=ylabel, output_folder=out_folder, V=V, parent_sim=parent_sim)

def closest_point(points,point):
    points = np.asarray(points)
    deltas = points - point
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)

def create_seed_coords(xmin, xmax, dx, ymin, ymax, dy, separation=7):
    """Creates n_seeds 2D coordinates on the grid defined by the bounds, with a minimum separation. e.g. to seed inital points for Monte Carlo."""
    coords=[]
    x=np.arange(xmin,xmax,dx)
    y=np.arange(ymin,ymax,dy)
    nx=len(x)
    ny=len(y)
    nx=int((nx%separation)/2)
    ny=int((ny%separation)/2)
    x=x[nx::separation]
    y=y[ny::separation]
    for xx in x:
        for yy in y:
            coords.append([xx,yy])
    return coords

def create_layer_coords(seed_coord, layer_index, xmin, xmax, dx, ymin, ymax, dy, separation=7):
    """Creates n_seeds 2D coordinates on the grid defined by the bounds, with a minimum separation. e.g. to seed inital points for Monte Carlo."""
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
    coords = np.array([coord for coord in coords if (xmin+dx<coord[0]<xmax-dx and ymin+dy<coord[1]<ymax-dy)])
    return coords

def create_seed(seed_index, seed_coord, INITIAL_FIELD, xlabel, ylabel, output_folder, V, parent_sim):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    temp_conf = f'''parent_sim='{parent_sim}'
INITIAL_FIELD='{INITIAL_FIELD}'
[{xlabel}, {ylabel}]={seed_coord}
V={V}
'''
    confname=f'{0:03d}_{seed_index:03d}'
    for coord in seed_coord:
        confname+=f'_{coord:.3f}'
    # confname+=INITIAL_FIELD
    confname = os.path.join(output_folder,confname)
    with open(confname+'.conf', 'w') as f:
        f.write(temp_conf)

def create_seeds(seed_coords, INITIAL_FIELDS, xlabel, ylabel, output_folder, V, parent_sim):


    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for INITIAL_FIELD in INITIAL_FIELDS:
        for seed_index,seed_coord in enumerate(seed_coords):

            create_seed(seed_index=seed_index, seed_coord=seed_coord, INITIAL_FIELD=INITIAL_FIELD, xlabel=xlabel, ylabel=ylabel, output_folder=output_folder, V=V, parent_sim=parent_sim)

def create_germ(layer_coord, dx, dy, seed_index, seed_coord, layer_index, xlabel, ylabel, output_folder, V, parent_sim):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    parent_filename=parent_label(seed_coord, seed_index, layer_coord, layer_index, dx, dy)
    temp_conf = f'''parent_sim='{parent_sim}'
parent_filename='{parent_filename}'
[{xlabel}, {ylabel}]={list(layer_coord)}
V={V}
'''
    confname=f'{layer_index:03d}_{seed_index:03d}'
    for coord in layer_coord:
        confname+=f'_{coord:.3f}'
    confname = os.path.join(output_folder,confname)
    with open(confname+'.conf', 'w') as f:
        f.write(temp_conf)

def create_germs(germ_coords, dx, dy, seed_index, seed_coord, layer_index, xlabel, ylabel, output_folder, V, parent_sim):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for layer_coord in germ_coords:
        create_germ(layer_coord=layer_coord, dx=dx, dy=dy, seed_index=seed_index, seed_coord=seed_coord, layer_index=layer_index, xlabel=xlabel, ylabel=ylabel, output_folder=output_folder, V=V, parent_sim=parent_sim)

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
