from lib_plt import *
from mpl_toolkits.mplot3d import Axes3D

if len(sys.argv) < 2:
	print('Please supply name of folder')
	sys.exit()

folder = sys.argv[1]
data_location=os.path.join('data',folder)

data_length = np.size(glob.glob(data_location+'\*.conf'))

def Data(i,k):
    global data, omegas, dos, r

    config_module = import_path(glob.glob(data_location+'\*.conf')[i]) 
    module_dict, to_import = import_all(config_module)
    globals().update({name: module_dict[name] for name in to_import})

    data = np.load(glob.glob(data_location+'\*.npz')[i], allow_pickle=True)
    dos = data['dos']
        
    dos = dos.sum((2,3))
    
    indices=FindIndicesOfArray(omegas, omegaMin, omegaMax)
    
    omegas=omegas[indices]
    dos=dos[:,:,indices]   
    r=[k,0]

    [x,y]=centre(r,n_x,n_y) #centred coords
    
    dos = dos[x,y,:]
    return()
#
def Sort_data(Data, data_length):
    '''Returns a list to call Data by the magnitude of mu, rather than using range() in a loop'''
    order=[]
    
    for i in range(data_length): 
        Data(i,0)
        order.append(V)
    
    temp = {val: key for key, val in enumerate(order)}
    res = list(map(temp.get, sorted(order)))
    return res
##########################################################
######################## Main ############################
##########################################################
def Main():
    res = Sort_data(Data, data_length)
    res=res[Vlim:-Vlim:Vstep]
    Z=[]
    for k in range(num_plts):
        for i in res:    
            Data(i,k)
            Z.append(dos)
    Z=np.array(Z)

    vmax=np.amax(Z)
    vmin=np.amin(Z)
    dv=vmax-vmin

    fig = plt.figure()
    title=r'$-\frac{1}{\pi}\Im\hat{G}^R(\omega+i\epsilon;\mathbf{r},\mathbf{r})$'
    # fig.suptitle(title)
    for k in range(num_plts):

        colors = [cm.Spectral(i) for i in np.linspace(0, 1, data_length)]
        cmap='Spectral'

        contourHeight=-1

        ax= fig.add_subplot(2,1,k+1,projection='3d')   

        y=[]
        Z=[]
        for i in res:    
            Data(i,k)
            y.append(V)
            Z.append(dos)
            
        Z=np.array(Z)
        x=np.real(omegas)
        
        X,Y = np.meshgrid(x,y)
        
        offset=10**(-6)


        
        lower_2D_projection = 8
        min=vmin-vmin*lower_2D_projection
        ZZ=offset*Z+min
        vmax2=offset*vmax+min+contourHeight
        vmin2=offset*vmin+min+contourHeight
        
        ax.contourf(X, Y, ZZ+contourHeight, 50, cmap=cmap, antialiased=True,
        vmax=vmax2, vmin=vmin2)
        ax.contourf(X, Y, ZZ+contourHeight, 50, cmap=cmap, antialiased=True,
        vmax=vmax2, vmin=vmin2)
        ax.plot_surface(X, Y, Z, cmap=cmap, rstride=1, cstride=1, antialiased=True,
        vmax=vmax, vmin=vmin, alpha=1, linewidth=0)
        ax.plot_surface(X, Y, Z, cmap=cmap, rstride=1, cstride=1, antialiased=True,
        vmax=vmax, vmin=vmin, alpha=1, linewidth=0)
        ax.plot_surface(X, Y, Z, cmap=cmap, rstride=1, cstride=1, antialiased=True,
        vmax=vmax, vmin=vmin, alpha=1, linewidth=0)
        ax.plot_surface(X, Y, Z, cmap=cmap, rstride=1, cstride=1, antialiased=True,
        vmax=vmax, vmin=vmin, alpha=1, linewidth=0)
        # ax.set(title=(title+r'$|_{{\mathbf{{r}}={}}}$'.format(r)))
        ax.text2D(0.68, 0.82, title+r'$|_{{\mathbf{{r}}={}}}$'.format(r), transform=ax.transAxes)

        ax.set(xlabel=r'$\omega/t$')
        ax.set(ylabel=r'$V/t$')
        # ax.set(zlabel='LDOS')
        # ax.tick_params(left=False, labelleft=False)
        handles, labels = ax.get_legend_handles_labels()
        ncol = int(np.ceil(len(handles)/2))
        zmin,zmax=0, 0.5
        ax.set_zlim(zmin,zmax)

        ax.set_zlim(vmin+contourHeight, vmax)
        label_loc=[zmin,zmin+zmax/2,zmax]
        ax.set(zticks=label_loc)
        
        ax.w_xaxis.set_pane_color((1,1,1,1)) #white background
        ax.w_yaxis.set_pane_color((1,1,1,1)) 
        ax.w_zaxis.set_pane_color((1,1,1,1)) 
        
        # ax.view_init(20,-120)
        
    fig.set_size_inches(w=latex_width, h=1.5*latex_width) 

    plt.tight_layout()
    # fig.subplots_adjust(hspace=0.25)
    
    print('Done')
    return fig
###############################################
#################### Main #####################
###############################################
latex_width=4.7747
omegaMin,omegaMax=-4,10
Vlim,Vstep=20,1
num_plts=2
fig=Main()
plt.show()
fig.savefig('out/'+'Spectrum_3D_'+folder+'.pdf', bbox_inches = "tight")