from lib_plt import *

if len(sys.argv) < 2:
	print('Please supply name of configuration file')
	sys.exit()
    
data_location='data'
data_location=os.path.join(data_location,'superconducting_state')

configuration = sys.argv[1]
confname = configuration.split('.conf')[0]
config_module = import_path(os.path.join(data_location,configuration))    
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

def Data():
    global dos, omegas

    data = np.load(glob.glob(os.path.join(data_location,'DOS_'+confname+'.npz'))[0], allow_pickle=True)
    dos = data['dos']
    
    indices=FindIndicesOfArray(omegas, omega_min, omega_max)
    dos = dos[:,:,:,:,indices]
    omegas=omegas[indices]
    
    dos=dos.sum((2,3))
    return
##########################################################
########################## Main ##########################
##########################################################
def Main():
    fig = plt.figure()
    ax = fig.add_subplot(111)
        
    color_index=[4,0]
    i=0
    cycle=10
    colors = [cm.nipy_spectral(i) for i in np.linspace(0, 1, cycle)]
    rr=pt[0]
    r=centre(rr,n_x,n_y)
    x,y=np.real(omegas),dos[r[0],r[1],:]
    ax.plot(x,y, linestyle='solid', color=colors[color_index[i]],label=pt[i])
    ax.fill_between(x,y,0,color=colors[color_index[i]], alpha=0.4)#, clip_on=False)

    i=i+1
    rr=pt[1]
    r=centre(rr,n_x,n_y)
    x,y=np.real(omegas),dos[r[0],r[1],:]
    ax.plot(x,y, linestyle='solid', color=colors[color_index[i]],label=pt[i])
    for i in range(4,-1,-1):        
        
        ax.fill_between(x,y,0,color=colors[i], alpha=0.4)

    ax.set(xlabel=r'$\omega$')  
    ax.set(ylabel=r'Density of states')

    text_DOS = ('Spectrum'
    '\n'
    f'$\epsilon={np.imag(omegas[0])}$'
    )    

    fig.text(.22, .81, text_DOS,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center')
                           
    legend = fig.legend(loc="upper left", 
    fancybox=True, shadow=True, #prop=fontP,
    bbox_to_anchor=(0.14,0.58,1,0.2))
    
    fig.set_size_inches(w=latex_width, h=4.5) 

    return fig, ax

pt=[[5,0],[0,0]]
omega_min,omega_max=np.min(omegas),np.max(omegas)#-1.6,1.6
latex_width=4.7747

Data()
fig, ax = Main()
plt.show()
fig.savefig('out/'+'Spectrum_'+confname+'.pdf', bbox_inches = "tight")