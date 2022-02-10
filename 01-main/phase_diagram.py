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

def create_config_files(INITIAL_FIELDS, xlabel, ylabel, xmin, xmax, dx, ymin, ymax, dy, n_layers, separation, FILENAME, V):
    """Creates configuration files for the phase diagram"""

    out_folder = os.path.join(ROOT_DIR,CONF,FILENAME,OUT)
    silent=sys.argv[-1]
    ClearOutFolder(out_folder,silent=silent)

    seed_coords=create_seed_coords(xmin=xmin, xmax=xmax, dx=dx, ymin=ymin, ymax=ymax, dy=dy, separation=separation)

    create_seeds(seed_coords=seed_coords, INITIAL_FIELDS=INITIAL_FIELDS, xlabel=xlabel, ylabel=ylabel, output_folder=out_folder, V=V)

    tmp=int(separation/2)
    if n_layers<tmp:
        raise ValueError("n_layers must be greater than half the separation distance! Using n_layers=separation/2 instead.")
        n_layers=tmp

    for seed_index,seed_coord in enumerate(seed_coords):
        for layer_index in range(1,n_layers):
            germ_coords=create_layer_coords(seed_coord=seed_coord, layer_index=layer_index, xmin=xmin, xmax=xmax, dx=dx, ymin=ymin, ymax=ymax, dy=dy, separation=separation)
            create_germs(germ_coords=germ_coords, dx=dx, dy=dy, seed_index=seed_index, seed_coord=seed_coord, layer_index=layer_index, xlabel=xlabel, ylabel=ylabel, output_folder=out_folder, V=V)

