from lib import *

class PhaseDiagram():

    def __init__(self,model):

        """Model must be a function of various variables, which returns a (TightBinding/BogoliubovdeGennes) class"""

        self.model=model
        self.directory=DATA
        self.initial_name='initial_convergence'
        self.initial_title=''
        self.filename=''
        self.xlabel='Independent variable'
        self.ylabel='Dependent variable'
        
        self.friction = 0
        self.max_iterations = 100
        self.absolute_convergence_factor = 0.00001

        self.indep_vars=[]

        self.plotting_func=None

        self.png=False
    
    def set_png(self):
        self.png=True
    
    def set_directory_name(self,directory_name):
        FILENAME=sys.argv[0].split('.')[0]
        self.data=os.path.join(DATA,FILENAME,directory_name)
        self.fig=os.path.join(FIG,FILENAME,directory_name)
        for directory in [self.data,self.fig]:
            if not os.path.exists(directory):
                os.makedirs(directory)

    def __delete_files(self):
        files=glob.glob(os.path.join(self.data,'*.npz'))
        if len(files)>0:
            print('Data files found')
            if YesNo('Delete?'):
                for f in files:
                    os.remove(f)

    def set_xlabel(self,xlabel):
        self.xlabel=xlabel

    def set_ylabel(self,ylabel):
        self.ylabel=ylabel

    def set_friction(self,friction):
        self.friction=friction

    def set_max_iterations(self,max_iterations):
        self.max_iterations=max_iterations

    def set_absolute_convergence_factor(self,absolute_convergence_factor):
        self.absolute_convergence_factor=absolute_convergence_factor

    def set_dep_vars_func(self,func_returning_dep_vars):
        """Requires a func which takes an instance of BogoloiubovdeGennes or TightBinding and returns the dep variables of interest for the phase diagram."""
        self.func_returning_dep_vars=func_returning_dep_vars

    def set_indep_vars(self,indep_vars):
        self.indep_vars=indep_vars

    def set_plotting_func(self,plotting_func):
        """plotting_func must be a funciton of bdg which returns fig,ax"""
        self.plotting_func=plotting_func

    def save_data(self,bdg):
        self.filename=f'{self.i:04d}_{self.indep_vars[self.i]}'
        filename=os.path.join(self.data,self.filename+'.npz')
        with open(filename, 'wb') as f:
                cPickle.dump(bdg,f)

    def save_fig(self,fig):
        if self.png:
            filename=os.path.join(self.fig,self.filename+'.png')
        else:
            filename=os.path.join(self.fig,self.filename+'.pdf')
        plt.savefig(filename, bbox_inches = "tight", dpi=DPI)
        plt.close()

    def __iter__(self):

        return self

    def __next__(self):

        x = self.indep_vars[self.i]
        bdg = self.model(x)
        
        if self.i>0:
            bdg._hartree=self.bdg._hartree
            bdg._fock=self.bdg._fock
            bdg._gorkov=self.bdg._gorkov
            bdg._hubbard_indices=self.bdg._hubbard_indices
            bdg._anomalous_indices=self.bdg._anomalous_indices
            bdg.U_entries=self.bdg.U_entries
        
        bdg.self_consistent_calculation(friction=self.friction, max_iterations=self.max_iterations, absolute_convergence_factor=self.absolute_convergence_factor)

        self.y = self.func_returning_dep_vars(bdg)

        del bdg.eigenvalues
        del bdg.eigenvectors
        del bdg._hamiltonian
        del bdg._tb_ham
        del bdg._hubbard_u

        self.save_data(bdg)

        self.bdg = bdg

        self.i+=1

    def execute(self):

        self.__delete_files()

        iteration = iter(self)

        for self.i in tqdm.tqdm(range(len(self.indep_vars)),desc='Indep var', total=len(self.indep_vars),leave=False):

            next(iteration)

            if type(self.plotting_func)!=type(None):
                self.save_fig(
                        self.plotting_func(self.bdg)
                        )

    def plot(self):
        for self.i in tqdm.tqdm(range(len(self.indep_vars)),desc='Indep var', total=len(self.indep_vars),leave=False):

            self.filename=f'{self.i:04d}_{self.indep_vars[self.i]}'
            filename=os.path.join(self.data,self.filename+'.npz')
            bdg = np.load(filename, allow_pickle=True)

            if type(self.plotting_func)!=type(None):
                self.save_fig(
                        self.plotting_func(bdg)
                        )

class PlotMinimisedData():

    def __init__(self):

        """Model must be a function of various variables, which returns a (TightBinding/BogoliubovdeGennes) class"""

        self.directory=DATA
        self.initial_name='initial_convergence'
        self.initial_title=''
        self.filename=''
        self.xlabel='Independent variable'
        self.ylabel='Dependent variable'
        
        self.friction = 0
        self.max_iterations = 100
        self.absolute_convergence_factor = 0.00001

        self.indep_vars=[]

        self.plotting_func=None

        self.png=False
    
    def set_png(self):
        self.png=True
    
    def set_data_directories(self,directories):
        FILENAME=sys.argv[0].split('.')[0]
        temp=[]
        for directory in directories:
            temp.append(os.path.join(DATA,FILENAME,directory))
        self.directories=temp
    
    def set_directory_name(self,directory_name):
        FILENAME=sys.argv[0].split('.')[0]
        self.data=os.path.join(DATA,FILENAME,directory_name)
        self.fig=os.path.join(FIG,FILENAME,directory_name)
        for directory in [self.data,self.fig]:
            if not os.path.exists(directory):
                os.makedirs(directory)

    def set_xlabel(self,xlabel):
        self.xlabel=xlabel

    def set_ylabel(self,ylabel):
        self.ylabel=ylabel

    def set_plotting_func(self,plotting_func):
        """plotting_func must be a funciton of bdg which returns fig,ax"""
        self.plotting_func=plotting_func
    
    def minimised_data(self):
        names=[]
        variables=[]
        free_energies=[]
        for directory in self.directories:
            directory=glob.glob(os.path.join(directory,'*.npz'))
            for file in directory:
                bdg = np.load(file, allow_pickle=True)
                names.append(file)
                variable=file.split('/')[-1].split('.npz')[0].split('_')[1]
                variables.append(variable)
                free_energies.append(bdg.free_energy[-1])
        # Sort:

        # names=[x for _, x in sorted(zip(variables, names))]
        # free_energies=[x for _, x in sorted(zip(variables, free_energies))]
        # variables=[x for _, x in sorted(zip(variables, variables))]

        # Minimise free energy:
        names=np.array(names)
        unique_variables=np.array(list(set(variables)),dtype=float)
        unique_variables=[x for _, x in sorted(zip(unique_variables, unique_variables))]
        variables=np.array(variables,dtype=float)
        free_energies=np.array(free_energies)
        minimised_names=[]
        for variable in unique_variables:
            indices=np.where(variables==variable)[0]
            index=np.argmin(free_energies[indices])
            minimised_names.append(names[indices[index]])

        self.variables = unique_variables

        self.minimised_names=minimised_names

        self.save_data()
    
    def save_data(self):
        FILENAME=sys.argv[0].split('.')[0]
        for file in self.minimised_names:
            txt=file
            filename=file.split('/')[-1].split('.npz')[0]
            filename=filename.split('_')[-1]
            filename+='.txt'
            filename=os.path.join(self.data,filename)
            with open(filename, 'w') as f:
                f.write(txt)

    def save_fig(self,fig):
        if self.png:
            filename=os.path.join(self.fig,f'{self.i:04d}_{self.variables[self.i]}.png')
        else:
            filename=os.path.join(self.fig,f'{self.i:04d}_{self.variables[self.i]}.pdf')
        plt.savefig(filename, bbox_inches = "tight", dpi=DPI)
        plt.close()

    def plot(self):

        files=self.minimised_names
        n=len(files)
        for self.i,file in tqdm.tqdm(enumerate(files),desc='Indep var', total=len(self.indep_vars),leave=False):
        
            filename=f'{self.variables[self.i]}'
            
            filename=os.path.join(self.data,filename+'.txt')
            with open(filename) as f:
                filename = f.readlines()[0]
                bdg = np.load(filename, allow_pickle=True)

                if type(self.plotting_func)!=type(None):
                    self.save_fig(
                            self.plotting_func(bdg)
                            )

    
