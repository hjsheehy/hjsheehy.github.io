from lib_plt import *
import matplotlib.ticker as mtick

# if len(sys.argv) < 2:
	# print('Please supply name of folder')
	# sys.exit()

# folder = sys.argv[1]
folder = 'normal_state_varying_mu'
data_location=os.path.join('data',folder)

data_length = np.size(glob.glob(data_location+'\*.conf'))

def Data(i):
    global data, omega, dis_straight, dos_straight, dis_diagonal, dos_diagonal, dis_imp, exact_straight, x_axis

    config_module = import_path(glob.glob(data_location+'\*.conf')[i]) 
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})
    
    data = np.load(glob.glob(data_location+'\*.npz')[i], allow_pickle=True)
    dos, omega = data['dos'], data['omega']
    
    dos = dos.sum((2,3))

    dis=2*6
    r=impurity_locations[0]
    r=[-dis,0]
    [x0,y0]=centre(r,n_x,n_y) #centred coords
    r1=[dis,0]
    [x1,y0]=centre(r1,n_x,n_y)
    x_axis=list(range(r[0], r1[0]))


    line=Straight_line(y0,x0,x1) #line from centre to horizontal edge
    dis_straight=np.array(list(range(len(line)))) #distance from centre at each pt
    dos_straight=Linecut(dos, line)
    exact_line=[[x,0] for x in x_axis]
    # exact_line=[[x,0] for x in dis_straight]
    exact_straight=[Analytical_Friedel(r0, V, mu, t, omega) for r0 in exact_line]

    line=Diagonal_line(x0, x1) #line from centre to corner edge
    dis_diagonal=np.array(list(range(len(line)))) #distance from centre at each pt
    dos_diagonal=Linecut(dos, line)
    exact_line=[[x,x] for x in dis_diagonal]
    exact_diagonal=np.array([Analytical_Friedel(r0, V, mu, t, omega) for r0 in exact_line])

    return()
##########################################################
######################## Main ############################
##########################################################
def Main():
    
    fig, axs = plt.subplots(nrows=ceil(data_length/2),ncols=2)#, subplot_kw=dict(polar=True))
    
    i=0
    for ax_row in axs:
        for ax in ax_row:
            Data(i)
            i+=1
            x0 = dis_straight
            label_x = r'No. sites moving horizontally from $r_0$'
            label_y = (r'''Density of states
            ''')

            y0 = dos_straight
            y1 = exact_straight
            
            label_y0 = f'$N={n_x}$'
            label_y1 = 'Analytical'
                
            if i==1:
                # handles, labels = fig.get_legend_handles_labels()
                label_y1 = 'Analytical'
            else:
                label_y0 = None
                label_y1 = None

            ax1 = ax.twinx()

            color='tab:red'
            ix,iy=i%2,floor(i/2)
            ax.plot(x_axis, y0,color=color,marker='o',markersize=3,label=label_y0)
            ax.tick_params(axis='y', labelcolor=color)
            
            ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
            
            color = 'tab:blue'
            ax1.plot(x_axis, y1,color=color,marker='x',markersize=3,label=label_y1)  
            ax1.tick_params(axis='y', labelcolor=color)      

            plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)            

            txt = f'$\mu={mu:.2f}$'
            # plt.title(txt)
            
            # v_min = np.min(y0)
            # ax.text(-0.9, v_min, txt,
                     # {'bbox': dict(boxstyle="square", alpha=0.8, fc="white",
                                   # ec="none", pad=0.2)}, ha='right', va='bottom')
            v_max = np.max(y1)
            ax1.text(-1, v_max, txt,
                     {'bbox': dict(boxstyle="square", alpha=0.8, fc="white",
                                   ec="none", pad=0.2)}, ha='right', va='top')

    fig.set_size_inches(w=latex_width, h=2*latex_width)

    # handles, labels = ax1.get_legend_handles_labels()

    legend = fig.legend(loc="upper center", 
    #fancybox=True, shadow=True, #prop=fontP,
    bbox_to_anchor=(0.275,-0.0,0.05,0.2))
    
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel(label_x)
    # plt.ylabel(label_y)
    plt.title(label_y)

    # plt.tight_layout()
    return fig

fig=Main()
# plt.show()

fig.savefig('out/'+f'LDOS_line_Analytical_versus_numerical_mosaic'+'.pdf', bbox_inches = "tight")