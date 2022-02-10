from lib_plt import *
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

if len(sys.argv) < 2:
	print('Please supply name of folder')
	sys.exit()

folder = sys.argv[1]
data_location=os.path.join('data',folder)

# labels=np.array([
# ['''Normal_state_43_43_0_0.5_0.75_1.00_0_0_0_0_1.21''',r'Normal'],
# ['''Singlet_43_43_0_0.5_0.75_1.00_1_1_0_0_1.21''',r'Singlet'],
# ['''Unitary_triplet_43_43_0_0.5_0.75_1.00_1_1_0_0_1.21''',r'$\boldsymbol{\eta}=\mathbf{\hat{x}}$'],
# ['''Sudeep_INT_state_0.025''',r'Non-unitary triplet']])

labels=np.array([
['''Two_orbital_Singlet_43_43_0_0.05_0.075_1.00_1_0.35_1.21''',r'Multiorbital singlet'],
['''Non-unitary_triplet_43_43_0_0.05_0.075_1.00_1_-0.35_-1j_0_1.21''',r'Non-unitary triplet']])

data_length=len(labels)

columns =[]
n_col = 6

for i in range(data_length):
    confname=os.path.join(data_location,labels[i,0])
    confname_DOS=os.path.join(data_location,'DOS_'+labels[i,0])

    data = np.load(confname_DOS+'.npz', allow_pickle=True)
    dos=data['dos']
    
    dos = dos.sum((2,3))

    f = np.fft.fft2(dos, axes=(0,1), norm='ortho')
    f = np.fft.fftshift(f, axes=(0,1))
    
    abs_f = np.abs(f)
    rho_plus = np.real(f + np.flip(f, -1))
    rho_minus = np.real(f - np.flip(f, -1))
    
    columns.append([abs_f, rho_plus, rho_minus])
    
columns = np.array(columns)

config_module = import_path(confname+'.conf') 
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

epsilon=np.imag(omegas[0])
n_omega=len(omegas)

# max/min without central bright spot and lines:
r=centre([0,0],n_x,n_y)
r=(* r,)
vmax=[]
vmin=[]
for i in range(3):
    temp=np.copy(columns[:,i,:,:,int(n_omega/2):])
    temp[:,r[0],r[1],:]=temp[:,0,0,:]
    vmax.append(np.amax(temp))
    vmin.append(np.amin(temp))
print(vmax)
print(vmin)
# vmax=vmax_set
# vmin=vmin_set

def Data(omega_input, columns):
    index = FindNearestValueOfArray(np.real(omegas), omega_input)
    omega=omegas[index]
      
    row = columns[:,:,:,:,index]
    return omega, row
##########################################################
######################## Main ############################
##########################################################
def Main():
    n_rows = len(omegas_plt)
                
    xlabel='$k_x a$'
    ylabel='$k_y a$'
    
    title0 = r'Multiorbital singlet'
    title1 = r'Non-unitary triplet'
    subtitle0 = r'$|\text{FT}|$'
    subtitle1 = r'$\rho_+$'
    subtitle2 = r'$\rho_-$'
    txt1 = f'$\epsilon={epsilon:.2f}$'   
    
    cmp1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black","orange","white"])
    cmp2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black","blue","green","yellow","red"])
    cmp3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black","purple","yellow","cyan"])
   
    interpolation = 'none'
    extent = [-n_x/2,n_x/2,-n_y/2,n_y/2]

    fig, ax = plt.subplots(n_rows, n_col, sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    for i in range(0,len(omegas_plt)):
        omega_input = omegas_plt[i]
        omega, row = Data(omega_input, columns)
        eV=np.real(omega)
        txt0 = f'$\omega={eV:.2f}$'

        row[0,0,:,:]
        im0 = ax[i][0].imshow(
            row[0][0].T, extent=extent,
            vmin=vmin[0], vmax=vmax[0],
            interpolation=interpolation,
            cmap=cmp1)    
        ax[i][0].text(0.39,0.11, 
            txt0,
            {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                           ec="none", pad=0.2)}, ha='center', va='center', transform=ax[i][0].transAxes)
        im1 = ax[i][1].imshow(
            row[0][1].T, extent=extent,
            vmin=vmin[1], vmax=vmax[1],
            interpolation=interpolation,
            cmap=cmp2)
        im2 = ax[i][2].imshow(
            row[0][2].T, extent=extent,
            vmin=vmin[2], vmax=vmax[2],
            interpolation=interpolation,
            cmap=cmp3)
        im3 = ax[i][3].imshow(
            row[1][0].T, extent=extent,
            vmin=vmin[0], vmax=vmax[0],
            interpolation=interpolation,
            cmap=cmp1)
        im4 = ax[i][4].imshow(
            row[1][1].T, extent=extent,
            vmin=vmin[1], vmax=vmax[1],
            interpolation=interpolation,
            cmap=cmp2)
        im5 = ax[i][5].imshow(
            row[1][2].T, extent=extent,
            vmin=vmin[2], vmax=vmax[2],
            interpolation=interpolation,
            cmap=cmp3)
        for j in range(n_col):
            ax[i][j].set_xticklabels([])
            ax[i][j].set_yticklabels([])
        
    text_x, text_y = 0.5,1.2
    for j in range(2):
        ax[0][0+3*j].text(text_x+0.015, text_y,
            subtitle0, ha='center', va='center', transform=ax[0][0+3*j].transAxes)

        ax[0][1+3*j].text(text_x, text_y,
            subtitle1, ha='center', va='center', transform=ax[0][1+3*j].transAxes)

        ax[0][2+3*j].text(text_x, text_y, 
            subtitle2, ha='center', va='center', transform=ax[0][2+3*j].transAxes)
                               
    ax[0][1].text(text_x, text_y+0.25,
        title0, ha='center', va='center', transform=ax[0][1].transAxes)

    ax[0][1+3].text(text_x, text_y+0.25,
        title1, ha='center', va='center', transform=ax[0][1+3].transAxes)


    axins0 = inset_axes(ax[-1][3],
                    width="10%", # width = 10% of parent_bbox width
                    height="70%", # height : 50%
                    loc=6)
    cbar0=plt.colorbar(im3, cax=axins0, ticks=[vmin[0],vmax[0]])
    cbar0.ax.set_yticklabels(['Low','High'], color='white')

    axins1 = inset_axes(ax[-1][4],
                    width="10%", # width = 10% of parent_bbox width
                    height="70%", # height : 50%
                    loc=6)
    cbar1=plt.colorbar(im4, cax=axins1, ticks=[vmin[1],vmax[1]])
    cbar1.ax.set_yticklabels(['Low','High'], color='white')


    axins2 = inset_axes(ax[-1][5],
                    width="10%", # width = 10% of parent_bbox width
                    height="70%", # height : 50%
                    loc=6)
    cbar2=plt.colorbar(im5, cax=axins2, ticks=[vmin[2],vmax[2]])
    cbar2.ax.set_yticklabels(['Low','High'], color='white')

    plt.tick_params(left=False,
                    bottom=False,
                    labelleft=False,
                    labelbottom=False)
                
#############################################
    # ax3.set(xlabel=xlabel)
    # ax4.set(xlabel=xlabel)
    # ax1.set(ylabel=ylabel)
    # ax3.set(ylabel=ylabel)
    
    # plt.setp(ax2.get_yticklabels(), visible=False)
    # plt.setp(ax4.get_yticklabels(), visible=False)
    # plt.setp(ax1.get_xticklabels(), visible=False)
    # plt.setp(ax2.get_xticklabels(), visible=False)
    
    # ax3.set_xticks([-n_x/2,0,n_x/2])
    # ax3.set_xticklabels(['$-\pi$',0,'$\pi$'])
    # ax3.set_yticks([-n_y/2,0,n_y/2])
    # ax3.set_yticklabels(['$-\pi$',0,'$\pi$'])
    # ax4.set_xticks([-n_x/2,0,n_x/2])
    # ax4.set_xticklabels(['$-\pi$',0,'$\pi$'])
    # ax1.set_yticks([-n_y/2,0,n_y/2])
    # ax1.set_yticklabels(['$-\pi$',0,'$\pi$'])
    
    # ax1.text(0.16,0.87, 
        # text_DOS,
             # {'bbox': dict(boxstyle="square", alpha=0.6, fc="white",
                           # ec="none", pad=0.2)}, ha='center', va='center', transform=ax1.transAxes)

    fig.set_size_inches(w=latex_width, h=1.7*latex_width) 

    # plt.tight_layout()
    
    return fig

latex_width=1.2*4.7747
# vmax_set=0.05
# vmin_set=-vmax_set
dOmega=0.05

for ii in range(4):
    omega_min, omega_max=0.5*ii, 0.5*(ii+1)
    omegas_plt=np.arange(omega_min, omega_max, dOmega)
    fig=Main()
    # plt.show()
    fig.savefig('out/'+f'QPI_large_mosaic_'+folder+f'_({ii+1}).pdf', bbox_inches = "tight")
print('Done')