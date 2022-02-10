from lib_plt import *
from matplotlib.offsetbox import AnchoredText
# from itertools import cycle 

if len(sys.argv) < 2:
	print('Please supply name of folder')
	sys.exit()
    
folder = sys.argv[1]
data_location=os.path.join('data',folder)

data_length = np.size(glob.glob(data_location+'\*.conf'))

def Data(i):  
    global data, omegas, x0, y0, label      
    confname=glob.glob(data_location+'\*.conf')[i]
    
    config_module = import_path(glob.glob(data_location+'\*.conf')[i]) 
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})
    
    label=(confname.split(folder+'\\')[1]).split('.conf')[0]
    
    data = np.load(glob.glob(data_location+'\*.npz')[i], allow_pickle=True)
    dos = data['dos']
    
    dos = dos.sum(3) #sum over orbital component (which is 1 for all simulations being plotted)

    indices=FindIndicesOfArray(omegas, omega_min, omega_max)
    dos = dos[:,:,:,indices]
    omegas=omegas[indices]
    
    [x0,y0]=centre(r,n_x,n_y) #centred coords
   
    dos = dos[x0,y0]
        
    dos_total=np.sum(dos,0)
    dos_up=dos[0]
    dos_down=dos[1]
    
    data=[dos_up, dos_down,dos_total]
    label=[r'$\uparrow$',r'$\downarrow$','total']
    return
##########################################################
########################## Main ##########################
##########################################################
def Main():
    fig, axs = plt.subplots(data_length)
    
    xlabel=r'$\omega$'
    ylabel=r'Density of states'
    
    # linewidth = np.linspace(1, 4, data_length)
    lines = ["-","-",":"]
    colors = ['k','r','b']
    
    j=1
    for ax in axs:
        Data(j)
        j-=1
        for i in range(3):
            y=data[i]
            x=np.real(omegas)
            
            # ax.plot(x, y,
            ax.plot(x, y, 
            label=label[i],
            linestyle=lines[i], color=colors[i])
                                   
            text_sigma = (f'$\sigma={2*np.imag(omegas[0])}$')    
            axs[j].text(2.2, 0.01, text_sigma,
                     {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                                       ec="none", pad=0.2)}, ha='center', va='center')


    for ax in axs.flat:
        ax.set(xlabel=xlabel, ylabel=ylabel)
        ax.label_outer()
        
    text_DOS = (f'''Spectrum
    $\mathbf{{r}}={r}$''')    
    axs[0].text(.24, .75, text_DOS,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                               ec="none", pad=0.2)}, ha='center', va='center')

    legend = axs[0].legend(loc="upper left", #framealpha=0.8,
    fancybox=True, shadow=True, 
    bbox_to_anchor=(0.02,0.63,1,0.2))
    
    fig.set_size_inches(w=latex_width, h=1.3*latex_width) 
    
    plt.tight_layout()
    return fig, axs
 
r=[0,5] 
omega_min,omega_max=0,2.5
latex_width=4.7747

fig,axs=Main()
plt.show()

fig.savefig('out/'+'Spectrum_spin_resolved_'+folder+'.pdf', bbox_inches = "tight")