from plt_lib import *
########################################################
########################## Data ########################
########################################################
layer=0
omega=0
n_var=4
n_legend=3
n_plt_data=len(filenames)
plt_data=np.zeros([n_var,n_legend,n_plt_data])
x=np.zeros([n_legend,50])
for i in range(n_plt_data):
    filename = filenames[i]
    globals().update(conf_file(filename))

    if n_orbitals==1 and n_spins==1:
        j=0
    if n_orbitals==2 and n_spins==1:
        j=1
    if n_orbitals==2 and n_spins==2:
        j=2
    if n_orbitals==1 and n_spins==2:
        continue

    data = np.load(filename, allow_pickle=True)
    dos, exec_time, mem = data['dos'], data['exec_time'], data['mem']
    ldos = LDOS(dos, omegas, omega, trace_over=True, layer=layer)
    n=n_x
    ldos=np.sum(ldos)/n**2
    exec_time=exec_time/3600
    mem=max(mem)*1048576 # MiB to b
    mem=mem*(10**(-9)) # b to gb
    k=int(n/2)-1
    plt_data[:,j,k]=[n,ldos,exec_time,mem]
########################################################
########################## Plot ########################
########################################################
def main():
    marker_list=['.','x','+']
    size_list=[5,7,6]

    n=int(np.max(plt_data[0,:,:]))
    n_list = np.arange(1,n+2,2)
    analytical = Analytical_Friedel([n,0], 0, mu, t, omega)
    analytical_y = np.ones(np.shape(n_list))*analytical
    x0 = plt_data[0]
    label_x = r'$n$'

    y0 = plt_data[1]
    label_y0 = r'$\text{{LDOS}}\,(r_n)$'
    
    y2 = plt_data[2]
    y3 = plt_data[3]
    label_y2= r'Execution time (h)'
    label_y3= r'Memory (gb)'      
        
    fig = plt.figure()
    ax0 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex= ax0)    
    ax3 = ax2.twinx()

    for i in range(n_legend):
        color='black'
        ax0.scatter(x0[i], y0[i],color=color,marker=marker_list[i],s=3*size_list[i],label=label_y0)
        color='tab:orange'
        ax2.scatter(x0[i], y2[i],color=color,marker=marker_list[i],s=3*size_list[i],label=label_y2)
        color = 'tab:blue'
        ax3.scatter(x0[i], y3[i],color=color,marker=marker_list[i],s=3*size_list[i],label=label_y3)  
    
        color='black'
        ax0.plot(n_list, analytical_y,color=color, linestyle='solid',linewidth=0.5)
        color='tab:blue'
        n_spins=1
        n_orbitals=1
        if i==1:
            n_orbitals=2
            n_spins=1
        if i==2:
            n_orbitals=2
            n_spins=2
        mem=[]
        for n in n_list:
            mem.append(memory([n,n,1,n_spins,n_orbitals])*(10**(-9))
)
        ax3.plot(n_list, mem,color=color, linestyle='solid',linewidth=0.5)

    color='black'
    ax0.set_ylabel(label_y0, color=color)
    ax0.tick_params(axis='y', labelcolor=color)
    
    color='tab:orange'
    ax2.set_ylabel(label_y2,color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_xlabel(label_x)
    
    color = 'tab:blue'
    ax3.tick_params(axis='y', labelcolor=color)
    ax3.set_ylabel(label_y3, color=color)
    # ax0.set_ylim([0.025,0.125])
    ax2.set_yscale("log")
    ax3.set_yscale("log")
    nmax=np.amax(plt_data[0,:])
    listOf_Xticks = np.arange(3, nmax+1, 4)
    plt.xticks(listOf_Xticks)
    
    plt.setp(ax0.get_xticklabels(), visible=False)
    
    labels = ['Memory without overhead', 'Spinless', 'Multiorbital spinless', 'Multorbital with spin']
    
    color='black'

    analytical = mlines.Line2D([], [], color=color, marker='None', linestyle='solid',
                          markersize=size_list[0], label=labels[0])
    blue_star = mlines.Line2D([], [], color=color, marker=marker_list[0], linestyle='None',
                          markersize=size_list[0], label=labels[1])
    red_square = mlines.Line2D([], [], color=color, marker=marker_list[1], linestyle='None',
                          markersize=size_list[1], label=labels[2])
    purple_triangle = mlines.Line2D([], [], color=color, marker=marker_list[2], linestyle='None',
                          markersize=size_list[2], label=labels[3])

    ax2.legend(handles=[analytical, blue_star, red_square, purple_triangle])

    plt.tight_layout()

    fig.set_size_inches(w=LATEX_WIDTH, h=4.5) 
    return fig
########################################################
######################### Caption ######################
########################################################
def caption():
    fermi_vector = Fermi_vector(mu, t, omega)
    friedel_wavelength = Friedel_wavelength(fermi_vector)

    text=rf'''Benchmarking to observe attenuation of finite size effects with
increasing n.
'''
    with open(output+'.txt', 'w') as f:
        f.write(text)

##########################################################
########################## Main ##########################
########################################################## 
show_then_save(main, caption)
