from plt_lib import *
########################################################
########################## Data ########################
########################################################
n_data=len(filenames)
n_list=[]
iterations_list=[]
exec_time_list=[]
hartree_list=[]
gorkov_list=[]
for i in range(n_data):
    ################ import ################
    filename = filenames[i]
    globals().update(conf_file(filename))

    data = np.load(filename, allow_pickle=True)
    hartree, gorkov, exec_time, mem, iterations = data['hartree'], data['gorkov'], data['exec_time'], data['mem'], data['iterations']
    n_list.append(n_x)
    hartree_list.append(hartree[:,-1])
    gorkov_list.append(gorkov[:,-1])
    iterations_list.append(iterations)
    exec_time_list.append(exec_time)

########################################################
########################## Plot ########################
########################################################
def main():
    label_x = r'$n$'
    label_y1 = r'$\Delta_\text{singlet}$'
    label_y2 = r'$\phi_\sigma$'
    x0=np.array(n_list)
    y0=np.transpose(gorkov_list)[0]
    y11=np.transpose(hartree_list)[0]
    y12=np.transpose(hartree_list)[1]
    y2=np.array(iterations_list)
    y3=np.array(exec_time_list)*(1/60)
    
    inds = x0.argsort()
    x0=x0[inds]
    y0=y0[inds]
    y11=y11[inds]
    y12=y12[inds]
    y2=y2[inds]
    y3=y3[inds]

    i=0
    label_x = r'$n_x=n_y$'

    label_y0 = r'$\Delta_{{\uparrow,\downarrow}}/t$'

    label_y1= r'$\phi_{{\uparrow,\downarrow}}/t$'

    label_y2= r'Iterations'    
    
    label_y3= r'Execution time (min)'    
    
    fig = plt.figure()

    color='tab:red'
    ax0 = fig.add_subplot(211)
    ax0.plot(x0, y0,color=color,marker='.',markersize=4,label=label_y0)
    ax0.set_ylabel(label_y0, color=color)
    ax0.tick_params(axis='y', labelcolor=color)

    ax1 = ax0.twinx()

    color = 'tab:blue'
    ax1.plot(x0, y11,color=color,marker='x',markersize=4,label=label_y1)  
    ax1.plot(x0, y12,color=color,marker='+',markersize=4,label=label_y1)  
    ax1.set_ylabel(label_y1, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    color='tab:orange'
    ax2 = fig.add_subplot(212, sharex= ax0)    
    ax2.plot(x0, y2,color=color,marker='.',markersize=4,label=label_y2)
    ax2.set_ylabel(label_y2,color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_xlabel(label_x)
    
    ax3 = ax2.twinx()

    color = 'tab:green'
    ax3.plot(x0, y3,color=color,marker='x',markersize=4,label=label_y3)  
    ax3.tick_params(axis='y', labelcolor=color)
    ax3.set_ylabel(label_y3, color=color)

    nmax=np.amax(n_list)
    listOf_Xticks = np.arange(3, nmax+1, 4)
    plt.xticks(listOf_Xticks)
    
    plt.setp(ax0.get_xticklabels(), visible=False)

    plt.tight_layout()

    fig.set_size_inches(w=latex_width, h=4.5) 

    return fig
########################################################
######################### Caption ######################
########################################################
def caption():
    #fermi_vector = Fermi_vector(mu, t, omega)
    #friedel_wavelength = Friedel_wavelength(fermi_vector)

    text=rf'''
'''
    with open(output+'.txt', 'w') as f:
        f.write(text)

##########################################################
########################## Main ##########################
########################################################## 
main()
plt.savefig(output+'.pdf', bbox_inches = "tight")
caption()
