from lib_plt import *

# if len(sys.argv) < 2:
	# print('Please supply name of configuration file')
	# sys.exit()

data_location = '../Data/Normal/'

conf_nsc_list=[]
conf_sc_list=[]
conf_nsc_list.append('spinless_V=1.21_varying_n/Normal_varying_n_79_-3.57_0_0_1.00_0_1.21e+00.conf')
conf_sc_list.append('self-consistent_spinless_V=1.21_n=79_mu=varying/Self-consistent_singlet_79_-3.57_0_0_1.00_1.21e+00_1.20e+00_0.conf')
conf_nsc_list.append('spinless_V=1.21_n=79_varying_mu/Normal_state_varying_V_79_-3.82_0_0_1.00_0_1.21e+00.conf')
conf_sc_list.append('self-consistent_spinless_V=1.21_n=79_mu=varying/Self-consistent_singlet_79_-3.82_0_0_1.00_1.21e+00_1.20e+00_0.conf')
conf_nsc_list.append('spinless_V=1.21_n=79_varying_mu/Normal_state_varying_V_79_-3.96_0_0_1.00_0_1.21e+00.conf')
conf_sc_list.append('self-consistent_spinless_V=1.21_n=79_mu=varying/Self-consistent_singlet_79_-3.96_0_0_1.00_1.21e+00_1.20e+00_0.conf')
conf_nsc_list.append('Normal_state_varying_V_79_-3.99_0_0_1.00_0_1.21e+00.conf')
conf_sc_list.append('self-consistent_spinless_V=1.21_n=79_mu=varying/Self-consistent_singlet_79_-3.99_0_0_1.00_1.21e+00_1.20e+00_0.conf')

def Data(configuration):
    global x_axis, dis_straight, dos_straight, dis_diagonal, dos_diagonal, exact_straight, exact_diagonal
    
    confname = configuration.split('.conf')[0]
    config_module = import_path(os.path.join(data_location,configuration))    
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})

    data = np.load(os.path.join(data_location,confname+'.npz'), allow_pickle=True)
    dos = data['dos']
    
    dos = dos.sum((2,3,4,5))

    index = FindNearestValueOfArray(np.real(omegas), omega)

    dos = dos[:,:,index]

    r=impurity_locations[0]
    [x0,y0]=centre(r,n_x,n_y) #centred coords
    x0=0
    x_axis=list(range(int(-n_x/2), int(n_x/2+1)))

    line=Straight_line(y0,x0,n_x) #line from centre to horizontal edge
    dis_straight=np.array(list(range(len(line)))) #distance from centre at each pt
    dos_straight=Linecut(dos, line)
    exact_line=[[x,0] for x in x_axis]
    # exact_line=[[x,0] for x in dis_straight]
    # exact_straight=[Analytical_Friedel(r0, V, mu, t, omega) for r0 in exact_line]
    exact_straight=[Analytical_LDOS(r0, V, mu, t, omega+1.0j*epsilon) for r0 in exact_line]
    
    line=Diagonal_line(x0, n_x) #line from centre to corner edge
    dis_diagonal=np.array(list(range(len(line)))) #distance from centre at each pt
    dos_diagonal=Linecut(dos, line)
    exact_line=[[x,x] for x in dis_diagonal]
    exact_diagonal=np.array([Analytical_Friedel(r0, V, mu, t, omega) for r0 in exact_line])
    return
##########################################################
######################## Main ############################
##########################################################
def Main():
    fig = plt.figure()
    ax0 = fig.add_subplot(111)
    # arrow_loc = []
    for j in range(len(conf_nsc_list)):
        configuration = conf_nsc_list[j]
        Data(configuration)
        
        x0 = x_axis
        label_x = r'No. sites horizontally away from impurity'
        label_y = r'Local density of states'

        y0 = dos_straight
        label_y0 = f'Non-self-consistent numerical model'
        
        y1 = exact_straight
        label_y1 = 'First-order perturbation analytically'

        arrow_loc = (x0[-26], (y0[-1]+y1[-1])/2)
        txt_loc = (x0[-14], y0[-1]*(1-0.24))
        print(txt_loc)
        ax0.annotate(f'$\mu={mu:.2f}$', xy=arrow_loc, xytext=txt_loc,
            arrowprops=dict(arrowstyle="->"))
        
        configuration = conf_sc_list[j]
        Data(configuration)
        
        y2 = dos_straight
        label_y2 = f'Self-consistent numerical model'

        color='tab:red'
        ax0.plot(x_axis, y0,color=color,marker='o',markersize=4)
        color='tab:green'
        ax0.plot(x_axis, y2,color=color,marker='.',markersize=4)
        ax0.set_xlabel(label_x)
        ax0.set_ylabel(label_y)
        color = 'tab:blue'
        ax0.plot(x_axis, y1,color=color,marker='x',markersize=4)  

    color='tab:red'
    ax0.plot(x_axis, y0,color=color,marker='o',markersize=4,label=label_y0)
    color='tab:green'
    ax0.plot(x_axis, y2,color=color,marker='.',markersize=4,label=label_y2)
    ax0.set_xlabel(label_x)
    ax0.set_ylabel(label_y)
    # ax0.tick_params(axis='y', labelcolor=color)
    av=np.average(y1[int(n_x/2)+2:])
    max=np.max(y1[int(n_x/2)+2:])
    max=max*1.002
    height=max-av
    min=av-height
    # min=0.98*np.min(y0[int(n_x/2)+2:])
    # ax0.set_ylim([min,max]) 
    
    color = 'tab:blue'
    ax0.plot(x_axis, y1,color=color,marker='x',markersize=4,label=label_y1)  

    fig.set_size_inches(w=latex_width, h=4.5)

    fig.legend(bbox_to_anchor=(0.82,1.16), fancybox=True, shadow=True)

    plt.tight_layout()
    return fig
    

omega = 0.0

fig=Main()
plt.show()
fig.savefig('out/'+f'LDOS_line_Analytical_versus_numerical'+'.pdf', bbox_inches = "tight")