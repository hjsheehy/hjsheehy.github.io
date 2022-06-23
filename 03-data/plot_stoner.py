import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import glob

directories=['Stoner_theory_Delta','Stoner_theory_Us']

names=[]
for directory in directories:
    names+=glob.glob(os.path.join(directory,'*'+'.npz'))
names_label=np.array([str(name.split('/')[-1].split('_')[-1]).split('.npz')[0] for name in names])
unique_labels=list(set(names_label))
names=[]
for directory in directories:
    for label in unique_labels:
        names+=glob.glob(os.path.join(directory,'*_'+label+'*.npz'))
    names_x=np.array([float(name.split('/')[-1].split('_')[0]) for name in names])
    names_y=np.array([str(name.split('_')[-2].split('.')[0]) for name in names])
    names=[x for _, x in sorted(zip(names_x, names))]
    names_y=[x for _, x in sorted(zip(names_x, names_y))]
    names_x=[x for _, x in sorted(zip(names_x, names_x))]
    names_i=np.arange(len(names))
    unique_x=list(set(names_x))

min_names=[]
for x in unique_x:
    indices=names_i[names_x==x]
    names_at_x=np.array(names)[indices]
    free_energy=[]
    for name in names_at_x:
        [x,y,z] = np.load(name, allow_pickle=True)
        x=np.array(x)
        y=np.array(y)
        x=x[:len(y[:,0])]
        free_energy.append(y[:,-1][-1])
    index=np.argmin(free_energy)
    min_name=np.array(names_at_x)[index]
    min_names.append(min_name)
names=min_names

#####################################
# temp
#####################################

fig, ax = plt.subplots(1,1,sharex='col')
field_index=0
minus_field_index=None
plus_field_index=None
absolute=False

n=len(np.unique(names_x))
markers=['>','<','.']
color = cm.gist_rainbow(np.linspace(0, 1, n))
xx=[]
yy=[]
zz=[]
fr=[]
free_energy=[]

s=2
markevery=1
for j,name in enumerate(names):
    i=list(set(names_x)).index(names_x[j])
    [x,y,z] = np.load(name, allow_pickle=True)
    x=np.array(x)
    y=np.array(y)
    x=x[:len(y[:,0])]
    if type(minus_field_index)!=type(None):
        y=y[:,field_index]-y[:,minus_field_index]
    elif type(plus_field_index)!=type(None):
        y=y[:,field_index]+y[:,plus_field_index]
    else:
        y=y[:,field_index]

    if absolute:
        y=np.abs(y)

    if names_y[j]=='min':
        c = color[i]
        marker=markers[2]
    elif names_y[j]=='fwd':
        # c = lighten_color(color[i],amount=0.3)
        c = color[i]
        marker=markers[0]
        ax.plot(x,y,marker=marker,markersize=s,color=c,label=f'${z:.2f}$',markevery=markevery)
    elif names_y[j]=='rev':
        # c = lighten_color(color[i],amount=0.3)
        c = color[i]
        marker=markers[1]
        ax.plot(x,y,marker=marker,markersize=s,color=c,markevery=markevery)

plt.show()
exit()
