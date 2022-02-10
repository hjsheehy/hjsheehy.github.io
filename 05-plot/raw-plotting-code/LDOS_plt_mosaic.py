from lib_plt import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

data_location=os.path.join('data','normal_state')

def Data():
    global dos_list, index
    config_module = import_path(os.path.join(data_location,confname[0]+'.conf'))    
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})

    dos_list=[]
    for i in range(4):
        data = np.load(glob.glob(os.path.join(data_location,'DOS_'+confname[i])+'.npz')[0], allow_pickle=True)
        dos = data['dos']
        dos = dos.sum((2,3))
        index = FindNearestValueOfArray(np.real(omegas), omega)
        dos = dos[:,:,index]
        dos_list.append(dos)
    return
##########################################################
########################## Main ##########################
##########################################################
def Main():
    eV=np.real(omegas[index])
    epsilon=np.imag(omegas[index])
    
    text1 = (r'a)')    
    text2 = (r'b)')
    text3 = (r'c)')
    text4 = (r'd)')

    xlabel='$x/a$'
    ylabel='$y/a$'
    # xticks=yticks=['$-\pi/2$',0,'$\pi/2$']
    
    xticks=[-int(n_x/2),0,int(n_x/2)]
    yticks=[-int(n_y/2),0,int(n_y/2)]
    xticks_labels=[f'${-int(n_x/2)}$',0,f'${int(n_x/2)}$']
    yticks_labels=[f'${-int(n_y/2)}$',0,f'${int(n_y/2)}$']
    fig = plt.figure()
    ax1 = fig.add_subplot(221)

    text_DOS = ('DOS map'
    '\n'
    f'$\omega={eV:.2f}$'
    '\n'
    f'$\epsilon={epsilon}$'
    )    
    
    vmax=np.amax(dos_list)
    vmin=np.amin(dos_list)


    interpolation = 'none'
    #
    # vmin, vmax = np.amin(dos_map), np.amax(dos_map)
    # dv=vmax-vmin
    # vmax=vmin+dv#/50
    #
    extent = [-n_x/2,n_x/2,-n_y/2,n_y/2]
    #
    i=0
    dos_map=dos_list[i]
    im=ax1.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)    
    # ax1.text(0.81,0.05, 
        # text1,
             # {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           # ec="none", pad=0.2)}, ha='center', va='center', transform=ax1.transAxes)
    
    i+=1
    dos_map=dos_list[i]
    ax2 = fig.add_subplot(222, sharey=ax1)
    ax2.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)
    # ax2.text(0.81,0.05, 
        # text2,
             # {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           # ec="none", pad=0.2)}, ha='center', va='center', transform=ax2.transAxes)

    i+=1
    dos_map=dos_list[i]
    ax3 = fig.add_subplot(223, sharex=ax1, sharey=ax1) 
    ax3.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)
    axins = inset_axes(ax3,
                    width="5%", # width = 10% of parent_bbox width
                    height="80%", # height : 50%
                    loc=6)
    cbar=plt.colorbar(im, cax=axins, ticks=[vmin,vmax])
    vmin,vmax=round(vmin,2),round(vmax,2)
    # cbar.ax.set_yticklabels([vmin,vmax], color='white')
    cbar.ax.set_yticklabels(['Low','High'], color='white')
    # ax3.text(0.79,0.05, 
        # text3,
             # {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           # ec="none", pad=0.2)}, ha='center', va='center', transform=ax3.transAxes)

    i+=1
    dos_map=dos_list[i]
    ax4 = fig.add_subplot(224, sharex=ax2, sharey=ax1)
    im=ax4.imshow(
        dos_map.T, extent=extent,
        vmin=vmin, vmax=vmax,
        interpolation=interpolation,
        cmap=newcmp)
    # ax4.text(0.77,0.05, 
        # text4,
             # {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           # ec="none", pad=0.2)}, ha='center', va='center', transform=ax4.transAxes)
    
    arrow_length=4
    arrow_alpha=1#0.8
    head_length=1
    head_width=1
    ax1.arrow(0.1,0+arrow_length+1,0,-arrow_length,
    fc="k", ec="k", alpha=arrow_alpha,
    head_width=head_width, head_length=head_length )
    
    ax1.arrow(5.1,0+arrow_length+1,0,-arrow_length,
    fc="g", ec="g", alpha=arrow_alpha,
    head_width=head_width, head_length=head_length )

    ax3.set(xlabel=xlabel)
    ax4.set(xlabel=xlabel)
    ax1.set(ylabel=ylabel)
    ax3.set(ylabel=ylabel)
    
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    
    ax3.set_xticks(xticks)
    ax3.set_xticklabels(xticks_labels)
    ax3.set_yticks(yticks)
    ax3.set_yticklabels(yticks_labels)
    ax4.set_xticks(xticks)
    ax4.set_xticklabels(xticks_labels)
    ax1.set_yticks(yticks)
    ax1.set_yticklabels(yticks_labels)
    
    ax1.text(0.19,0.86, 
        text_DOS,
             {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax1.transAxes)

    fig.set_size_inches(w=latex_width, h=latex_width) 

    plt.tight_layout()
    
    print('Done')
    return fig
confname=[
'Normal_state_43_43_-3.57_0_0_1.00_0_1.21e+00',
'Normal_state_two_impurities_43_43_-3.57_0_0_1.00_0_1.21e+00',
'Normal_state_corral_43_43_-3.57_0_0_1.00_0_1.21e+00',
'Normal_state_random_43_43_-3.57_0_0_1.00_0_1.21e+00']

omega=0
Data()

newcmp=cm.YlOrRd
omega=0
fig = Main()
plt.show()
fig.savefig('out/'+'LDOS_mosaic_normal_state'+'.pdf', bbox_inches = "tight")