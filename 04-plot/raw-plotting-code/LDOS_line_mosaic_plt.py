from lib_plt import *
import matplotlib.ticker as mtick

if len(sys.argv) < 2:
	print('Please supply name of folder')
	sys.exit()

folder = sys.argv[1]
# data_location=os.path.join('data',folder)
data_location=os.path.join('../Data/Normal/',folder)

data_length = np.size(glob.glob(data_location+'\*.conf'))

def Data(i):
    global dos, omegas, index, dis_straight, dos_straight, dis_diagonal, dos_diagonal, dis_imp

    config_module = import_path(glob.glob(data_location+'\*.conf')[i]) 
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})
    
    data = np.load(glob.glob(data_location+'\*.npz')[i], allow_pickle=True)
    dos = data['dos']
    
    dos = dos.sum((2,3,4,5))
    
    index = FindNearestValueOfArray(np.real(omegas), omega)

    dos = dos[:,:,index]
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
##########################################################
######################## Main ############################
##########################################################
def Main():
    label_sites = r'No. sites from $r_0$ to $r_\text{horiz}$'
    label_diagonal = r'No. sites from $r_0$ to $r_\text{diag}$'
    label_impurity_straight = r'Distance $r_0$ to $r_\text{horiz}$'
    label_impurity_diagonal = r'Distance $r_0$ to $r_\text{diag}$'
    
    colors = [cm.rainbow(i) for i in np.linspace(0, 1, data_length)]

    fig = plt.figure()

    ax1 = fig.add_subplot(221)
    for i in range(data_length): 
        Data(i)
        ax1.plot(dis_straight, dos_straight, label=f'${n_x}$', color=colors[i])
        
    ax2 = fig.add_subplot(222, sharey=ax1)
    for i in range(data_length):
        Data(i)
        ax2.plot(dis_imp, dos_straight, label=f'${n_x}$', color=colors[i])

    ax3 = fig.add_subplot(223, sharex=ax1, sharey=ax1)
    for i in range(data_length):
        Data(i)
        ax3.plot(dis_diagonal, dos_diagonal, label=f'${n_x}$', color=colors[i])
    handles, labels = ax3.get_legend_handles_labels()

    ax4 = fig.add_subplot(224, sharex=ax2, sharey=ax1)
    for i in range(data_length):
        Data(i)
        ax4.plot(dis_imp, dos_diagonal, label=f'${n_x}$', color=colors[i])
        
    lower_limit = 0 
    upper_limit = np.max([np.max(dos_straight), np.max(dos_diagonal)])*3.1
    # plt.ylim(lower_limit, upper_limit)

    ax1.set(xlabel=label_sites)
    ax2.set(xlabel=label_impurity_straight)
    ax3.set(xlabel=label_diagonal)
    ax4.set(xlabel=label_impurity_diagonal)
    # ax1.set(ylabel='LDOS')
    # ax3.set(ylabel='LDOS')
    fig.text(-0.01, 0.5, 'Local density of states', va='center', rotation='vertical')
    
    ax2.xaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax4.xaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)
    ncol = int(np.ceil(len(handles)/2))
    # ncol = len(handles)
    # fontP = FontProperties()
    # fontP.set_size('xx-small')
    legend = fig.legend(handles, labels, title='Number of sites in each direction, $n\equiv n_x=n_y$', loc="upper center", mode = "expand", ncol = ncol, 
    fancybox=True, shadow=True, #prop=fontP,
    bbox_to_anchor=(0,0.95,1,0.2))

    # k_F = Fermi_vector(mu, t, omega)
    # friedel_wavelength = Friedel_wavelength(k_F)
    # text_fermi = (f'Friedel_wavelength $={friedel_wavelength:.2f}$')
    # fig.text(0.81, .14, text_fermi,
             # {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           # ec="none", pad=0.2)}, ha='center', va='center')


    # plt.subplots_adjust(wspace=-.5, hspace=-1.)

    fig.set_size_inches(w=latex_width, h=5) 

    plt.tight_layout()
    
    print('Done')
    return fig
    
omega = 0.0
latex_width=4.7747
fig=Main()
Data(0)
print(f'Omega={np.real(omegas[index]):.3f}')
plt.show()
fig.savefig('out/'+'LDOS_line_mosaic_'+folder+'.pdf', bbox_inches = "tight")
# fig.savefig('LDOS_line_mosaic_'+folder+'.pdf', bbox_inches = "tight")