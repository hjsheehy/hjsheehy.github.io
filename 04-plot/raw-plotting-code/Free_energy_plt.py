from lib_plt import *
from scipy.sparse.linalg import eigsh
# from matplotlib.ticker import MaxNLocator

if len(sys.argv) < 2:
	print('Please supply name of folder')
	sys.exit()
    
folder = sys.argv[1]
data_location=os.path.join('data',folder)

data_length = np.size(glob.glob(data_location+'\*.conf'))

i=0
confname=glob.glob(data_location+'\*.conf')[i]

config_module = import_path(glob.glob(data_location+'\*.conf')[i]) 
module_dict, to_import = import_all(config_module)
globals().update({name: module_dict[name] for name in to_import})

def Data(i, Delta, phiUp, phiDown):  
    global w
    # Delta=dDelta*i
    n_x=n_y=n
    n_sites = n_x * n_y
    dof = n_sites * n_spins * n_orbitals
    
    onsite_spin=-mu*np.eye(n_spins)-phiUp*np.diag([1,0])-phiDown*np.diag([0,1])
    onsite_orbital=np.eye(n_orbitals)

    SC_spin=1.0j*Pauli_y
    SC_orbital=-Delta*np.eye(n_orbitals)   
    
    h0 = H0(onsite_spin, onsite_orbital, nn_spin, nn_orbital, im_spin, im_orbital,
            impurity_locations, n_x, n_y)
    mf = Mean_field(SC_spin, SC_orbital, n_sites)
    ham = H(h0, mf)
    w,v = la.eigh(ham)
    correlator=On_site_correlator(v,n_x,n_y,dof,n_sites,n_spins,n_orbitals)
    # Xi=np.sum(correlator[:,:,1,0,:,:dof]-correlator[:,:,0,1,:,:dof])/(2*n_orbitals)
    Xi=np.sum(correlator[:,:,1,0,:,:dof])/(n_orbitals)
    # density=Density_anomalous(v, dof)
    # dm=Density_matrix(density,n_x,n_y,dof,n_sites,n_spins,n_orbitals)
    Xi=np.sum(np.abs(correlator))/(2*n_orbitals)
    # print(Xi)
    
    density=Density_normal(v, dof)
    dm=Density_matrix(density,n_x,n_y,dof,n_sites,n_spins,n_orbitals)
    NUp=np.sum(dm[:,:,1,:,:dof])/(n_orbitals)
    NDown=np.sum(dm[:,:,0,:,:dof])/(n_orbitals)
    # print(Xi)
    
    # print(np.sum(dm[0,0,0,0,dof:]))
    
    # onsite_spin=-mu*np.eye(n_spins)
    # w,v = la.eigh(h0)
    # density=Density_normal(v, dof)
    # dm=Density_matrix(density,n_x,n_y,dof,n_sites,n_spins,n_orbitals)
    # E0=np.sum(dm[:,:,:,:,:dof])/(n_sites*n_spins*n_orbitals)

    return NUp, NDown, Xi
##########################################################
########################## Main ##########################
##########################################################
def old_Main():
    global x, y1, y2, y3, y4
    x=np.zeros(data_length)
    y1=np.zeros(data_length)
    y2=np.zeros(data_length)
    y3=np.zeros(data_length)
    y4=np.zeros(data_length)
    y1[0]=phiUpInitial
    y2[0]=phiDownInitial
    y3[0]=DeltaInitial
    
    Data(0,phiUpInitial,phiDownInitial,DeltaInitial)
    E_g=(-mu-2*t)#-np.sum(np.where(w < 0, w, 0))/dof       
    print(E_g)
    
    y4[0]=E_g
    for i in range(1,data_length):
        # energy=np.sum(w[:dof])
        # y[i]=energy
        # x[i]=Delta
        phiUp, phiDown, Delta = Data(i, y1[i-1], y2[i-1], y3[i-1])
        # freeEnergy=-mu*(phiUp+phiDown)-U*(phiUp*phiDown+Delta*np.conj(Delta))
        E_g=(-mu-phiDown-2*t)-np.sum(np.where(w < 0, w, 0))/dof       
        freeEnergy=E_g        
        # freeEnergy=-U*(phiUp*phiDown+Delta*np.conj(Delta))+2*U*(phiUp*phiDown)+2*U*Delta*np.conj(Delta)+E_g
        y1[i] = (1-f)*phiUp+f*y1[i-1]
        y2[i] = (1-f)*phiDown+f*y2[i-1]
        y3[i] = (1-f)*Delta+f*y3[i-1]
        y4[i] = freeEnergy
        x[i]=i
    # print(y1[-1])
    # print(y2[-1])
    indices=x.argsort()
    x=x[indices]
    y1=y1[indices]
    y2=y2[indices]
    y3=y3[indices]
    y4=y4[indices]
    
def Main():
    global x, y1, y2, y3, y4
    phiUp = phiDown = 0
    x=np.zeros(data_length)
    y4=np.zeros(data_length)
    DeltaList = np.linspace(DeltaStart,DeltaStop,data_length,endpoint=True)
    for i in range(data_length):
        Delta = DeltaList[i]
        x[i] = Delta
        # y1[i] = phiUp 
        # y2[i] = phiDown
        # y3[i] = Delta 
        NUp, NDown, Xi = Data(i, phiUp, phiDown, Delta)
        E_g = (-mu-phiDown)+np.sum(np.where(w < 0, w, 0))/dof
        # freeEnergy=E_g-(U/n_sites)*(Xi*np.conj(Xi)+NUp*NDown)+phiUp*NUp+phiDown*NDown+Delta*np.conj(Xi)+np.conj(Delta)*Xi
        # freeEnergy=-(U/n_sites)*(Xi*np.conj(Xi))+Delta*np.conj(Xi)+np.conj(Delta)*Xi
        # Z = np.prod(1+np.exp(-w/T))
        # freeEnergy = -T*np.log(Z)
        freeEnergy = -T*np.sum(np.log(1+np.exp(-w/T)))
        y4[i] = freeEnergy

T=0.5
data_length=35
DeltaStart,DeltaStop=-1,1
U=2.8
mu=-2.3
n=7
V=0
#f=0.8
Main()

ax = plt.figure().gca()
plt.plot(x,y4,color='red',marker='o')
# plt.xlim([0,data_length])
# plt.ylim([0,6.5])
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.show()