from lib_plt import *
# from itertools import cycle 

if len(sys.argv) < 2:
	print('Please supply name of folder')
	sys.exit()
    
folder = sys.argv[1]
data_location=os.path.join('..\\Data',folder)
data_length = np.size(glob.glob(data_location+'\*.conf'))
print(data_length)
# labels=np.array([
# ['''Normal_state_43_43_0_0.05_0.075_1.00_0_0_0_0_1.21''',r'Normal'],
# ['''Singlet_43_43_0_0.05_0.075_1.00_1_1_0_0_1.21''',r'Singlet'],
# ['''Two_orbital_Singlet_43_43_0_0.05_0.075_1.00_1_0.35_1.21''',r'Multiorbital singlet'],
# ['''Unitary_triplet_43_43_0_0.05_0.075_1.00_1_1_0_0_1.21''',r'$\boldsymbol{\eta}=\mathbf{\hat{x}}$'],
# ['''Unitary_triplet_43_43_0_0.05_0.075_1.00_1_0_1_0_1.21''',r'$\boldsymbol{\eta}=\mathbf{\hat{y}}$'],
# ['''Unitary_triplet_43_43_0_0.05_0.075_1.00_1_0_0_1_1.21''',r'$\boldsymbol{\eta}=\mathbf{\hat{z}}$'],
# ['''Non-unitary_triplet_43_43_0_0.05_0.075_1.00_1_-0.35_-1j_0_1.21''',r'Non-unitary triplet']])

labels=np.array([
['''Normal_state_43_43_0_0.05_0.075_1.00_0_0_0_0_1.21''',r'Normal'],
['''Singlet_43_43_0_0.05_0.075_1.00_1_1_0_0_1.21''',r'Singlet'],
['''Two_orbital_Singlet_43_43_0_0.05_0.075_1.00_1_0.35_1.21''',r'Multiorbital singlet'],
['''Non-unitary_triplet_43_43_0_0.05_0.075_1.00_1_-0.35_-1j_0_1.21''',r'Non-unitary triplet']])

# labels=np.array([
# ['''LaNiGa2_N''',r'Normal, V=1.21'],
# ['''LaNiGa2_SC_V=0''',r'Singlet, V=0'],
# ['''LaNiGa2_SC_V=1.21''',r'Singlet, V=1.21']])

data_length=len(labels)

def Data(i):  
    global dos, omegas, x0, y0
    confname=os.path.join(data_location,labels[i,0])
    confname_DOS=os.path.join(data_location,labels[i,0])
    
    config_module = import_path(confname+'.conf') 
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})
    
    # label=(confname.split(folder+'\\')[1]).split('.conf')[0]
    
    data = np.load(confname_DOS+'.npz', allow_pickle=True)
    dos = data['dos']
    
    # print(np.shape(dos))

    dos = dos.sum(2) #sum over spin component (which is 1 for all simulations being plotted)
    # dos = dos.sum(2) #sum over spin component (which is 1 for all simulations being plotted)
    
    indices=FindIndicesOfArray(omegas, omega_min, omega_max)
    dos = dos[:,:,:,indices]
    omegas=omegas[indices]
    
    [x0,y0]=centre(r,n_x,n_y) #centred coords
    return
##########################################################
########################## Main ##########################
##########################################################
def Main():
    fig = plt.figure()
    ax = fig.add_subplot(111)

    linewidth = [1,1,1,3,2,1,1]
    lines = ["-","-","-","--","-.",":","-"]
    colors = [cm.hsv(i) for i in np.linspace(0.1, 1, data_length)]
    
    # linecycler = cycle(lines)

    for i in range(data_length):
        Data(i) 
        
        # dos_total = dos[x0,y0,0,0,:] + dos[x0,y0,1,1,:]
        dos_total = dos[x0,y0,0,:] + dos[x0,y0,1,:]
        
        # dos_up = dos[x0,y0,0,:]
        # dos_down = dos[x0,y0,1,:]
        # y=dos[x0,y0,0,:] #Takes only one spin component (since singlet mean-field is degerenate anyway)
        x=np.real(omegas)
        
        # ax.plot(x, y,
        ax.plot(x, dos_total, 
        # label=r'${\boldsymbol\eta}=$'+f'${np.round(eta,2)}$, $s={round(s,2)}, \delta={round(delta,3)}$, $\Delta={round(Delta,2)}$', 
        label=labels[i,1],
        linestyle=lines[i], color=colors[i], linewidth=linewidth[i])
        # ax.plot(x, dos_up, linestyle=lines[1], color=colors[i])
        # ax.plot(x, dos_down, linestyle=lines[2], color=colors[i])

    # ax.set_ylim([0,0.6])

    ax.set(xlabel=r'$\omega$')  
    ax.set(ylabel=r'Density of states')

    text_DOS = (f'''Spectrum
    $\epsilon={np.imag(omegas[0])}$
    $\mathbf{{r}}={r}$''')    

    fig.text(.21, .81, text_DOS,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center')
                           
    fig.set_size_inches(w=latex_width, h=4.5) 

    legend = fig.legend(loc="upper center", mode = "expand", ncol = 3, 
    fancybox=True, shadow=True, #prop=fontP,
    bbox_to_anchor=(0,0.87,1,0.2))

    return fig, ax
    
latex_width=4.7747  
r=[5,0] 
omega_min,omega_max=-2.5,2.5
fig,ax=Main()
plt.show()

# fig.savefig('out\\'+'Spectrum_comparison_zero_chem_pot.pdf', bbox_inches = "tight")
fig.savefig('out\\'+'Spectrum_comparison_nonunitary_VS_multiorbital_0.1.pdf', bbox_inches = "tight")
# fig.savefig('out/'+'Spectrum_comparison_'+'nonunitary_VS_multiorbital'+'.pdf', bbox_inches = "tight")