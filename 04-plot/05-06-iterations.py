from plt_lib import *
########################################################
########################## Data ########################
########################################################
i=0
linestyles=['solid','dashed']
    
fig = plt.figure()
axs = fig.add_subplot(111)

for filename in filenames:
    globals().update(conf_file(filename))

    import _pickle as cPickle
    with open(filename, 'rb') as f:
        data = cPickle.load(f)

    energy=0

    centre=data.centre
    gorkov=data.gorkov_iterations
    DeltaUp=gorkov[0]
    DeltaDn=gorkov[1]

    title=f'Self-consistent gradient descent\n$U_T={U_T}, U_R={U_R}, \mu={mu}, V={V}$'

    y0=DeltaUp
    y1=DeltaDn

    linestyle=linestyles[i]
    plt.plot(y0,label=r'$\Delta_{\uparrow\uparrow}^{+-}(\mathbf{r})$'+f', {state}', linestyle=linestyle)
    plt.plot(y1,label=r'$\Delta_{\downarrow\downarrow}^{+-}(\mathbf{r})$'+f', {state}', linestyle=linestyle)

    label_x = r'Iterations'
    label_y = r'Gorkov mean fields'

    axs.set_xlabel(label_x)
    axs.set_ylabel(label_y)
    axs.set_title(title)

    axs.legend()

    i+=1

fig.set_size_inches(w=LATEX_WIDTH, h=5) 
plt.savefig(output+f'INT_Unitary_iterations.pdf', bbox_inches = "tight")
