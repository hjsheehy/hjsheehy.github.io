from lib_plt import *

if len(sys.argv) < 2:
	print('Please supply name of configuration file')
	sys.exit()
    
data_location='..\\Data'

configuration = sys.argv[1]
confname = configuration.split('.conf')[0]

config_module = import_path(os.path.join(data_location,configuration))    
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

def Data():
    global dos_line, distance, line, omegas, r0, r1

    data = np.load(os.path.join(data_location,confname+'.npz'), allow_pickle=True)
    dos = data['dos']
    
    indices=FindIndicesOfArray(omegas, omega_min, omega_max)
    dos = dos[:,:,:,:,:,:,indices]
    omegas=omegas[indices]
    
    dos=dos.sum((2,3))
    dos=dos[:,:,0,0]+dos[:,:,1,1]
    
    [x0,y0]=centre(r0,n_x,n_y) #centred coords
    [x1,y1]=centre(r1,n_x,n_y)
    line=Straight_line(y0,x0,x1) #line from centre to horizontal edge
   
    ### Transpose line:
    # line=np.array(line)[:,[1,0]]

    dis_straight=np.array(list(range(len(line)))) #distance from centre at each pt

    line=Diagonal_line(x0, x1) #line from centre to corner edge
    r0=list(np.array(line[0])-np.array([x0,y0]))
    r1=list(np.array(line[-1])-np.array([x0,y0]))
    dis_diagonal=np.sqrt(2) * np.array(list(range(len(line)))) #distance from centre at each pt
    dos_diagonal=np.array(Linecut(dos, line))

    # dis_imp = [2*x/n_x for x in range(len(line))]
    
    dos_line = dos_diagonal
    distance = dis_diagonal
    return
##########################################################
########################## Main ##########################
##########################################################
# def Linecut_plot(data, positions):
    # '''Density of states of the Bogoliubov quasiparticles as a function of 
    # energy.'''
    # fig, ax = plt.subplots(figsize=(8,6.5))
    # ax.plot(omega, dos_linecut, c='red', linestyle='solid', linewidth=2.0)
    # ax.set(xlabel=r'$\omega$', ylabel=r'LDOS(\omega)')

                           
    # plt.show()
    # return fig, ax
#
def Main():
    cycle=len(distance)
    # colors = [cm.nipy_spectral(i) for i in np.linspace(0, 1, cycle)]
    colors = [cm.plasma(i) for i in np.linspace(0, 1, cycle)]
    fig = plt.figure()
    ax = fig.add_subplot(111)

    max=np.max(dos_line)
    min=np.min(dos_line)
    norm=np.max([abs(min),abs(max)])
    norm=1
    for i in range(len(distance)-1,-1,-1): #length of line, inclusive slice
        vert=distance[i]+dos_line[i,:]/norm
        x=np.real(omegas)
        y=vert
        color=colors[i%cycle]
        ax.plot(x, y, color=color, clip_on=False)
        ax.fill_between(x,y,0,color=color, alpha=0.4, clip_on=False)
    ax.set(xlabel=r'$\omega$')  
    ax.set(ylabel=r'Distance along line $(a)$')
    ax.set(ylim=[distance[0]-0.5*1.1, distance[-1]+0.5*1.3])

    text_DOS = ('Linecut'
    '\n'
    f'Start: {r0}'
    '\n'
    f'End: {r1}'
    '\n'
    f'$\epsilon={np.imag(omegas[0])}$'
    )    
    
    ax.set_title(r'$-\frac{1}{\pi}\Im\hat{G}^R(\omega+i\epsilon;\mathbf{r},\mathbf{r})$', pad=25)

    fig.set_size_inches(w=latex_width, h=5.5) 

    fig.text(.23, .795, text_DOS,
             {'bbox': dict(boxstyle="square", alpha=0.8, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center')
    return fig, ax

latex_width=4.7747  

omega_min,omega_max=-5,5
r0=[0,0]
r1=[15,0]

Data()

fig, ax = Main()
plt.show()

fig.savefig('out\\Linecut_LaNiGa2_singlet.pdf', bbox_inches = "tight")