import numpy as np
#
def aniostropic_circle(Rx,Ry,N):
    phi = np.linspace(0,2*np.pi,N,endpoint=False)
    impurities=[]
    for i in range(N):
        impurities.append([Rx*np.cos(phi[i]),Ry*np.sin(phi[i])])
    return [np.rint(impurities).astype('int')]
##########################################################
########################## Main ##########################
##########################################################
Rx=Ry=14.5
N=15
x=aniostropic_circle(Rx,Ry,N)
print(x)