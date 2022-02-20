from lib_plt import *

if len(sys.argv) < 2:
	print('Please supply name of folder')
	sys.exit()

folder = sys.argv[1]
data_location=os.path.join('data',folder)

data_length = np.size(glob.glob(data_location+'\*.conf'))

def Data(i):
    global data, omegas, dis_straight, dos_straight, dis_diagonal, dos_diagonal, dis_imp

    config_module = import_path(glob.glob(data_location+'\*.conf')[i]) 
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})

    data = np.load(glob.glob(data_location+'\*.npz')[i], allow_pickle=True)
    dos = data['dos']
    
    index = FindNearestValueOfArray(np.real(omegas), omega)
    dos = dos[:,:,:,:,index]
    
    dos = dos.sum((2,3))
    
    # dos=dos/dof

    # r=impurity_locations[0]=[0,0]
    r=[0,0] #since no impurity in this model
    [x0,y0]=centre(r,n_x,n_y) #centred coords

    line=Straight_line(y0,x0,n_x) #line from centre to horizontal edge
    dis_straight=np.array(list(range(len(line)))) #distance from centre at each pt
    dos_straight=Linecut(dos, line)

    line=Diagonal_line(x0, n_x) #line from centre to corner edge
    dis_diagonal=np.sqrt(2) * np.array(list(range(len(line)))) #distance from centre at each pt
    dos_diagonal=Linecut(dos, line)
    
    dis_imp = [2*x/n_x for x in range(len(line))]

    return()
#
def Sort_data(Data, data_length):
    '''Returns a list to call Data by the magnitude of mu, rather than using range() in a loop'''
    order=[]
    
    for i in range(data_length): 
        Data(i)
        order.append(mu)
    
    temp = {val: key for key, val in enumerate(order)}
    res = list(map(temp.get, sorted(order)))
    return res
##########################################################
######################## Main ############################
##########################################################
def Main():
    colors = [cm.rainbow(i) for i in np.linspace(0, 1, data_length)]

    label_sites = r'Sites from impurity to centre-edge'

    fig = plt.figure()

    ax1 = fig.add_subplot(111)
    
    res = Sort_data(Data, data_length)
            
    j=0
    for i in res:    
        Data(i)
        ax1.plot(dis_straight, dos_straight, label=round(mu, 2), color=colors[j])
        j+=1
        
    ax1.set(xlabel=label_sites)
    ax1.set(ylabel='Local density of states')
    # ax1.tick_params(left=False, labelleft=False)
    handles, labels = ax1.get_legend_handles_labels()
    ncol = int(np.ceil(len(handles)/2))

    # sort both labels and handles by labels
    # labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: float(t[0])))
    
    legend = fig.legend(handles, labels, loc="lower right", ncol = 2,
    fancybox=True, framealpha=0.8, #prop=fontP,
    bbox_to_anchor=(0.45,0.15,0.5,0.5),
    title='Chemical potential $\mu$')

    # k_F = Fermi_vector(mu, t, omega)
    # Friedel_wavelength = Friedel_wavelength(k_F)
    # text_fermi = (f'Friedel_wavelength $={Friedel_wavelength:.2f}$')
    # fig.text(0.81, .14, text_fermi,
             # {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           # ec="none", pad=0.2)}, ha='center', va='center')
                           
    # text_av = (f'Mean value\nN=43: ${av[0]:.3e}$\nAnalytical: ${av[1]:.3e}$')
    # fig.text(0.76, .30, text_av,
             # {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           # ec="none", pad=0.2)}, ha='center', va='center')

    # plt.subplots_adjust(wspace=-.5, hspace=-1.)

    fig.set_size_inches(w=latex_width, h=4.5) 
    
    plt.tight_layout()
    
    print('Done')
    return fig
###############################################
#################### Main #####################
###############################################
latex_width=4.7747
omega=0
fig=Main()
plt.show()
fig.savefig('out/'+'LDOS_line_'+folder+'.pdf', bbox_inches = "tight")