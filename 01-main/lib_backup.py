#!/usr/bin/env python

# PythU Python self-consistent Hubbard U module
# January, 2022
__version__='1.0.0'
# 
# To add in next version: 
# -function of momentum for U, mean-fields and tight binding
# -function of space for placing impurities
# -fourier tranform these functional hamiltonians using fourier_transform_hamiltonian
# -add linecut to plot
# -change spin to spin_polarisation in plots (None, Up, Down)
# -implement 3D ARPES (Fermi surface)
# -automate name creation in conf.py
# -class Phase_diagram():
#       x-label
#       y-label
#       x_bound(min,max,delta)
#       y_bound(min,max,delta)
#       starting_pts.append('default',')
#       initial_mean_fields=a list of labels, e.g. ['int', 'unitary']
# -combine sweeps and choose minimum of free energy
# -automated method of using bulk mean-fields as an applied bias onto the surface layer
# -add dictionary to field records
############################################################################       
# Scanning tunnelling microscopy (STM), quasiparticle interference (QPI) and 
# angle-resolved photoemission spectrocopy (ARPES) simulation of (topological)
# non-unitary spin-triplet unconventional superconductivity and magnetism, in 
# quantum dots, thin films and 3D, with impurities.

# Copyright (C) 2022, Henry Joseph Sheehy
# 
# PythU is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PythU is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# A copy of the GNU General Public License should be available
# alongside this source in a file named TBC.  If not,
# see <http://www.gnu.org/licenses/>.
#
# PythU is availabe at http://www. TBC

import sys
import os
ROOT_DIR = os.path.dirname(
             os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT_DIR)
sys.path.append('src/')
from imports import *
###################################################################
############################ General ##############################
###################################################################
Pauli_z = np.array([[1,0],[0,-1]])
Pauli_y = np.array([[0,-1.0j],[1.0j,0]])
Pauli_x = np.array([[0,1],[1,0]])
Pauli_plus = (Pauli_x+1.0j*Pauli_y)*(1./2.)
Pauli_minus = (Pauli_x-1.0j*Pauli_y)*(1./2.)
Pauli_vec = np.dstack([Pauli_x, Pauli_y, Pauli_z])
Pauli_vec = np.moveaxis(Pauli_vec,-1,0)

ALL=slice(None)

def dagger(array):
    array = np.conj(array.T)
    return array

def wrap(positions, dimensions):
    temp = np.array(dimensions)/2
    return np.add(np.mod(np.add(positions, temp), dimensions), -temp)
###################################################################
############################ Indexing #############################
###################################################################
def Index(dimensions: tuple) -> "array":
    '''Returns an array of dimensions with elements
    given by the indexing Zn_Z'''
    return np.reshape(np.arange(np.prod(dimensions)),dimensions,'F')

def coordinates_to_indices(array, dimensions):
    dim=len(dimensions)
    temp=np.mod(array,dimensions)
    multipler=np.ones([dim],dtype=int)
    for i in range(1,dim):
        multipler[i:]*=dimensions[i-1]
    temp=np.multiply(temp,multipler)
    temp=np.sum(temp,axis=-1)
    return temp

###################################################################
##################### Analytical Friedel ##########################
###################################################################
def Fermi_vector(mu, t, omega):
    return np.sqrt((mu+omega)/t+4)
#
def Friedel_wavelength(k_F):
    return np.pi/k_F
#   
def Analytical_Friedel(r, V, mu, t, omega):
    m_eff = 1/(2*t)
    k_omega = Fermi_vector(mu, t, np.real(omega))
    rho = m_eff/(2*np.pi)
    const = V*(m_eff**2)/(2*np.pi)
    Bessel_J = j0(k_omega*la.norm(r))
    Bessel_Y = y0(k_omega*la.norm(r))
    pert = const*Bessel_J*Bessel_Y
    return rho + pert
#
def Free_Greens_Function(r, mu, t, omega):
    m_eff = 1/(2*t)
    k_omega = Fermi_vector(mu, t, omega)
    return -1.0j*(m_eff/2)*hankel1(0, k_omega*la.norm(r))
#
def Analytical_LDOS(r, V, mu, t, omega):
    zero = np.finfo(float).eps
    rho = Free_Greens_Function(zero, mu, t, omega)
    pert = V*(Free_Greens_Function(r, mu, t, omega))**2
    return -(1/np.pi)*np.imag(rho + pert)


def FindNearestValueOfArray(array, value):
    '''Returns index of nearest value in the array'''
    return (np.abs(array - value)).argmin()

def FindIndicesOfArray(array, bound1, bound2):
    '''Returns indices of array within boundary'''
    upperBound = max(bound1, bound2)
    lowerBound = min(bound1, bound2)
    return np.where(np.logical_and(array>=lowerBound, array<=upperBound))[0]
###################################################################
########################### Display array #########################
###################################################################
def DisplayArray(ax,array):
    array.real[abs(array.real)<1e-13]=0
    array.imag[abs(array.imag)<1e-13]=0
    array=np.real_if_close(array)

    ax.matshow(np.abs(array).T, cmap=plt.cm.Blues)
    [x,y]=np.shape(array)
    for i in range(x):
        for j in range(y):
            c = array[i,j]
            ax.text(i, j, str(c), va='center', ha='center')
    return ax
###################################################################
######################### Hamiltonian #############################
###################################################################
# class phase_diagram():
#     def __init__(self):
#         self.x_label='x'
#         self.y_label='y'
#         self.x_arange=[0, 5, 0.5]
#         self.y_arange=[0, 5, 0.5]
#         self.n_pts=[5]
# 
#     def set_mean_field():

class Orbital():
    """An atomic orbital, instantiated with a position and label"""

    def __init__(self, label):
        self.label=label

    def __repr__(self):
        return "Orbital()"

    def __str__(self):
        return "Orbital labelled {self.label}."

class Atom():
    """An atom, instantiated with a position, a label and an empty
    list of orbitals"""
    _orbitals=[]

    def __init__(self, position, label):
        self.position=position
        self.label=label

    def __repr__(self):
        return "Atom()"

    def __str__(self):
        
        return f"Atom lablled {self.label}, located at {self.position} with {self.n_orbitals} orbitals."

    def add_orbital(self, orbital):
        if not isinstance(orbital,Orbital):
            raise ValueError('orbital is not an instance of the Orbital class!')
        self._orbitals.append(orbital)

    @property
    def orbitals(self):
        return [orb.label for orb in self._orbitals]

    @property
    def n_orbitals(self):
        return len(self._orbitals)

def Crystal_lattice():
    """Crystal lattice"""
    def __init__(self, basis_vectors):
        basis_vectors=self.basis_vectors
    
    @property
    def n_dimensions(self):
        return len(basis_vectors)

class tight_binding():
    """Tight-binding model on a lattice defined by the tight_binding parent class
    
    Attributes
    ----------

    Parameters
    ----------
    conf : .conf file
        configuration file

    Attributes
    ----------
    dm : array
        normal density matrix
    adm : array
        anomalous density matrix
    t_dm : array
        normal thermal density matrix
    t_adm : array
        anomalous thermal density matrix
    """

    def __init__(self, dimensions, n_spins, basis=None, orbitals=None, pbc=None, k_space=False):
        
        self._model='tb'
        self.k_space = k_space

        self.dimensions=dimensions
        self.n_spins=n_spins
        self.basis=np.array(basis,dtype=float)
        self.orbitals=orbitals
        self.n_orbitals = len(orbitals)
        self.n_dim = len(dimensions)
        self.temperature = 0 

        self.hoppings=[]
        self.label_hoppings=[]
        self.impurities=[]
        
        self.inverse_basis = la.inv(self.basis)
        self.reciprocal_basis = 2*np.pi * self.inverse_basis

        self.hartree_iterations=None
        self.fock_iterations=None
        self.gorkov_iterations=None

        self.print_V=None
        self.print_V_mf=None
        self.print_Eg=None
        self.print_free_energy=None

        #self.k_dim = k_dim

        if type(pbc)==type(None):
            self.pbc=np.ones([self.n_dim],dtype=bool)
        else:
            self.pbc=pbc

        # avoid hopping to itself when one-dimensional
        # avoid double counting when two-dimensional
        for i,dimension in enumerate(dimensions):
            if dimension<=2:
                self.pbc[i]=False

        dim=self.dimensions
        self.edge = np.floor(0.5*np.array(dim),dtype=float)
        for i in range(len(dim)):
            if self.pbc[i]:
                self.edge[i]=+100000000

        dim=self.dimensions
        #self.edge = np.floor(0.5*np.array(dim),dtype=float)
        #for i in range(len(dim)):
        #    if self.pbc[i]:
        #        self.edge[i]=+100000000

        if len(basis)>self.n_dim:
            raise ValueError('Overcomplete basis!')
        for i, orb in enumerate(self.orbitals):
            if len(orb)<self.n_dim:
                raise ValueError(f'Orbital {orb} not {self.n_dim}D!')
        for i, vec in enumerate(basis):
            if len(vec)<self.n_dim:
                raise ValueError(f'Basis vector {vec} not {self.n_dim}D!')

        # Append a dimension for the orbitals:
        self.extended_dimensions  = np.append(self.dimensions, [self.n_spins,self.n_orbitals])
        self.n_cells = np.prod(self.dimensions)
        self.n_atoms = self.n_cells * self.n_orbitals
        self.n_dof = self.n_atoms * n_spins

        self.centre = np.array(np.array(self.dimensions)/2, dtype=int)  # centred-coordinates
        self.coord = np.reshape(np.indices(self.extended_dimensions),[len(self.extended_dimensions),self.n_dof],'F').T

        self.coord[:,:self.n_dim] = wrap(self.coord[:,:self.n_dim], self.dimensions)

        self.coord_cell = self.coord[:,:self.n_dim]

        self._position_cell = np.dot(self.coord_cell[:self.n_cells],basis)

        self.cartesian = np.array([[self._position_cell[i] + self.orbitals[m] for m in range(self.n_orbitals)] for i in range(self.n_cells)])
        self.reciprocal_cartesian = np.dot(self.cartesian,self.reciprocal_basis)

        self.coordinates_to_cell_position = np.reshape(self._position_cell, [*self.dimensions,self.n_dim],'F')
        self.coordinates_to_orbital_position = np.reshape(self.cartesian, [*self.dimensions,self.n_orbitals,self.n_dim],'F')
        self.reciprocal_to_orbital_position = np.reshape(self.reciprocal_cartesian, [*self.dimensions,self.n_orbitals,self.n_dim],'F')
        self.com=[]
        for i in range(self.n_dim):
            com=np.sum(np.array(self.orbitals)[:,i])/self.n_orbitals
            self.com.append(com)
        self.com=self.com
        ####################################################

        self._hamiltonian = None

    def set_temperature(self,temperature):
        self.temperature=temperature

    #################################### 
    ########### Hamiltonian ############
    ####################################

    def _onsite(self, onsite_amplitude, spin=None, orbital=None, axes=None, location=None):

        if type(orbital)!=type(None):
            if orbital>=self.n_orbitals:
                raise ValueError(f'Orbital index {orbital} greater than n_orbitals={self.n_orbitals}!')

        if np.isscalar(onsite_amplitude):
            if type(orbital)==type(None):
                orbitals = np.ones(self.n_orbitals)
            else:
                orbitals=np.zeros(self.n_orbitals)
                orbitals[orbital]=1
            orbit_tensor=np.diag(orbitals)

            if type(spin)==type(None):
                spins = np.ones(self.n_spins)
            else:
                spins=np.zeros([self.n_spins])
                spins[spin]=1
            spin_tensor=np.diag(spins)
            temp = onsite_amplitude*kron(spin_tensor,orbit_tensor)

        elif np.shape(onsite_amplitude)[0]==self.n_spins:
            if np.shape(onsite_amplitude)==np.shape(['spinUp','spinDown']):
                spin_tensor=np.diag(onsite_amplitude)
            elif np.shape(onsite_amplitude)==np.shape([['UpUp','UpDown'],['DownUp','DownDown']]):
                spin_tensor = onsite_amplitude
            else:   
                raise ValueError('onsite_amplitude must be scalar, pair or matrix (spins), or spin-orbit matrix')

            if type(orbital)==type(None):
                orbitals = np.ones(self.n_orbitals)
            else:
                orbitals=np.zeros([self.n_orbitals])
                orbitals[orbital]=1
            orbit_tensor=np.diag(orbitals)
            
            temp = kron(spin_tensor,orbit_tensor)

        elif np.shape(onsite_amplitude)[0]==self.n_spins*self.n_orbitals:
            temp = onsite_amplitude

        sites = np.zeros([self.n_cells])
        if type(location)==type(None):
            if type(axes)==type(None):
                axes=[ALL for i in range(self.n_dim)]
            axes=tuple(axes)
            index = Index(self.dimensions)
            index = (index[axes])
            sites[index] = 1

        else:
            sites[location]=1

        sites=np.diag(sites)

        return kron(temp, sites)

    def set_onsite(self, onsite_amplitude, spin=None, orbital=None, axes=None, location=None):
        if type(self._hamiltonian)==type(None):
            self._hamiltonian = np.zeros([self.n_dof, self.n_dof], dtype=complex_dtype) 
        self._hamiltonian += self._onsite(onsite_amplitude, spin, orbital, axes, location)
        
    def set_hopping(self, hopping_amplitude, orb_i=None, orb_f=None, hop_vector=None, label='', add_time_reversal=True):

        if type(self._hamiltonian)==type(None):
            self._hamiltonian = np.zeros([self.n_dof, self.n_dof], dtype=complex_dtype) 
        self._hamiltonian += self._hopping(hopping_amplitude=hopping_amplitude, orb_i=orb_i, orb_f=orb_f, hop_vector=hop_vector, add_time_reversal=add_time_reversal)
        
        if type(orb_i)==type(None):
            orb_i=0
        if type(orb_f)==type(None):
            orb_f=0

        self.hoppings.append([orb_i,orb_f,hop_vector])
        self.label_hoppings.append(label)

    def _hopping(self, hopping_amplitude, orb_i=None, orb_f=None, hop_vector=None, add_time_reversal=True):
        
        for orb in [orb_i,orb_f]:
            if type(orb)!=type(None):
                if orb>=self.n_orbitals:
                    raise ValueError(f'Orbital index {orb} greater than n_orbitals={self.n_orbitals}!')
            
        if type(hop_vector)==type(None):
            hop_vector=np.zeros([self.n_dim])

        if np.array_equal(hop_vector,np.zeros([self.n_dim])):
            temp = np.eye(self.n_cells)
        else:
            temp=np.zeros([self.n_cells,self.n_cells])
            x=np.indices(self.dimensions)
            x=np.moveaxis(x,0,-1)
            x=np.add(np.mod(np.add(x, self.centre), self.dimensions), -self.centre)
            y=np.copy(x)
            y=y+hop_vector
            # cut out hop outside of boundary:
            indices=np.invert(np.any(y>self.edge,axis=-1))
            x = x[indices]
            y = y[indices]
            x=coordinates_to_indices(x,self.dimensions)
            y=coordinates_to_indices(y,self.dimensions)
            x=x.flatten()
            y=y.flatten()
            temp[x,y]+=1
        
        if np.isscalar(hopping_amplitude):
            hopping_amplitude = hopping_amplitude*np.eye(self.n_spins)

        if np.shape(hopping_amplitude)[0]==self.n_spins:
            if (np.shape(hopping_amplitude)==np.shape(np.ones([self.n_spins,self.n_spins]))):
                spin_tensor = hopping_amplitude
            else:
                spin_tensor = np.diag(hopping_amplitude)

            if self.n_orbitals>1:
                if type(orb_i)==type(None) or type(orb_f)==type(None):
                    raise ValueError('orb_i, orb_f or spin_orbit_tensor not given!')
                else:
                    orbit_tensor=np.zeros([self.n_orbitals, self.n_orbitals])
                    orbit_tensor[orb_i,orb_f]=1
                    spin_orbit_tensor=kron(spin_tensor,orbit_tensor)
            else:
                spin_orbit_tensor = spin_tensor
        elif np.shape(hopping_amplitude)[0]==self.n_spins*self.n_orbitals:
            spin_orbit_tensor = hopping_amplitude

        temp=kron(spin_orbit_tensor,temp)

        onsite = bool(orb_i==orb_f and np.array_equal(hop_vector, np.zeros([self.n_dim])))
        if add_time_reversal and not onsite:
            temp+=dagger(temp)

        # print('spin tensor:')
        # print(spin_tensor)
        # print('orbital tensor:')
        # print(orbit_tensor)

        return temp

    def set_impurities(self, impurity_amplitude, locations, spins=None, orbitals=None):
        
        for location in locations:
            self.set_onsite(onsite_amplitude=impurity_amplitude, spins=spins, orbitals=orbitals, location=location)
        
        if type(spins)==type(None):
            spins=np.arange(self.n_spins)
        if type(orbitals)==type(None):
            orbitals=np.arange(self.n_orbitals)
        self.impurities.append([locations,list(spins),list(orbitals)])

    def set_magnetic_impurities(self, M, locations, orbitals=None):
        spin_matrix=np.einsum('i,ijk->jk',M,Pauli_vec)
        self.set_impurities(spin_matrix, locations, orbitals=orbitals)

    def unravel_hamiltonian(self):
        self._hamiltonian = np.reshape(self._hamiltonian, np.append(self.extended_dimensions,self.extended_dimensions), 'F')

    def ravel_hamiltonian(self):
        self._hamiltonian = np.reshape(self._hamiltonian, [self.n_dof,self.n_dof], 'F')
    def fourier_coeff(self):
        x=self.coord_cell
        k=np.copy(x)
        k=k*2*np.pi/(self.dimensions)-np.pi
        coeff=np.dot(x,k.T)
        coeff=np.exp(-1.0j*coeff)
        return coeff

    def fourier_transform_psi(self,v):
        coeff = self.fourier_coeff()
        if np.shape(v)[0]==2*self.n_dof:
            coeff=np.block([[coeff,np.conj(coeff)],[np.conj(coeff),coeff]])
        return np.dot(np.conj(coeff), v)

    def inv_fourier_transform_psi(self,v):
        coeff=np.conj(self.fourier_coeff())
        return np.dot(coeff, v)

    def fourier_transform_hamiltonian(self, transform=[], inverse_transform=[]):

        ravelled=False
        if np.shape(self._hamiltonian)==(self.n_dof,self.n_dof):
            ravelled=True
            self.unravel_hamiltonian()
        
        self._hamiltonian = np.fft.ifftn(self._hamiltonian,axes=transform)
        transform = tuple(np.array(transform)+len(self.extended_dimensions))
        self._hamiltonian = np.fft.fftn(self._hamiltonian,axes=transform)
    
        self._hamiltonian = np.fft.fftn(self._hamiltonian,axes=inverse_transform)
        inverse_transform = tuple(np.array(inverse_transform)+len(self.extended_dimensions))
        self._hamiltonian = np.fft.ifftn(self._hamiltonian,axes=inverse_transform)

        if ravelled:
            self.ravel_hamiltonian()
    
        ########################################
        ######### Statistical Mechanics ########
        ########################################

    def _density_matrix(self, eigenvectors):
        return np.einsum('in,in->in',eigenvectors[:self.n_dof],np.conj(eigenvectors[:self.n_dof]),optimize=True)
    
    def solve(self):
        t = time.time()
        w,v = la.eigh(self._hamiltonian, overwrite_a=True)
        #from scipy.sparse.linalg import eigsh
        # w,v = eigsh(ham, k=50, which='SM', return_eigenvectors=True)
        if self.k_space:
            v=self.inv_fourier_transform_psi(v)
            self.k_space=False

        self.exec_time = time.time() - t
        return w,v

    def plt_energy(self,fig,ax,axis,density_of_states):
        axes = np.arange(self.n_dim)
        axes = tuple(np.delete(axes,axis))
        y=np.sum(density_of_states,axes)
        y=np.fft.fftshift(y,axes=axis)
        interpolation = 'none'
        vmin, vmax = np.amin(y), np.amax(y)
        yy=np.real(self.energy_interval)
        ymin,ymax=min(yy),max(yy)
        x=np.shape(y)[0]
        extent = [-x/2,x/2,ymin,ymax]
        cmap = ListedColormap(cm.afmhot(np.linspace(0, 0.75, 256)))
        im = ax.imshow(y.T, interpolation=interpolation, aspect='auto', extent=extent, origin='lower',
        vmin=vmin, vmax=vmax,
        cmap=cmap)
        return fig,ax

class bogoliubov_de_gennes(tight_binding):
    """Bogoliubov-de Gennes model from TB parent class  
    
    Attributes
    ----------

    Parameters
    ----------
    conf : .conf file
        configuration file

    Attributes
    ----------
    dm : array
        normal density matrix
    adm : array
        anomalous density matrix
    t_dm : array
        normal thermal density matrix
    t_adm : array
        anomalous thermal density matrix
    """

    def __init__(self, dimensions, n_spins, basis=None, orbitals=None, pbc=None): #, k_dim=None):

        self.tb_model = tight_binding(dimensions, n_spins, basis, orbitals, pbc)
    
        super().__init__(dimensions, n_spins, basis, orbitals, pbc)

        self._model='bdg'

        self.iterations=0

        self.friction = 0.7
        self.max_iterations = 1000
        self.absolute_convergence_factor = 0.00001
        
        self.hartree=np.zeros([self.n_dof], dtype=complex_dtype)
        self.fock=np.zeros([self.n_dof,self.n_dof], dtype=complex_dtype)
        self.gorkov=np.zeros([self.n_dof,self.n_dof], dtype=complex_dtype)

        self._hartree_indices = []
        self._fock_indices = []
        self._gorkov_indices = []

        self._hubbard_u = None
        
        self._hartree_print_indices = []
        self._fock_print_indices = []
        self._gorkov_print_indices = []

        self.hartree_iterations=[]
        self.fock_iterations=[]
        self.gorkov_iterations=[]
        
        self.V=[]
        self.V_mf=[]
        self.Eg=[]
        self.free_energy=[]
        
        self.print_V=False
        self.print_V_mf=False
        self.print_Eg=False
        self.print_free_energy=False

    def set_mean_field_hamiltonian(self):

        # initialise tight binding hamiltonian if doesn't exist:
        try: 
            self._tb_ham
        except:
            self._tb_ham = self._hamiltonian

        n_dof = self.n_dof

        try:
            self._hubbard_indices
        except:
            print('Hubbard U tensor not given!')
            self._hubbard_u = np.diag(self.hartree)+self.fock+self.gorkov
            self._set_hubbard_indices()
            self.fock=self.fock[self._hubbard_indices]
            self.gorkov=self.gorkov[self._hubbard_indices]

        hartree=self.hartree
        fock=self.fock
        gorkov=self.gorkov

        hamiltonian = np.zeros([2*n_dof,2*n_dof],dtype=complex_dtype)
        hamiltonian[:n_dof,:n_dof]=self._tb_ham-np.diag(hartree)
        hamiltonian[n_dof:,n_dof:]=-(np.conj(self._tb_ham)-np.conj(np.diag(hartree)))
        hamiltonian[self._hubbard_indices[0],self.anomalous_indices[1]]=-fock
        hamiltonian[self.anomalous_indices[0],self._hubbard_indices[1]]=-(-np.conj(fock))
        hamiltonian[self._hubbard_indices[0],self.anomalous_indices[1]]=-np.conj(-gorkov)
        hamiltonian[self.anomalous_indices[0],self._hubbard_indices[1]]=-gorkov
        self._hamiltonian = hamiltonian

        #fig, ax = plt.subplots(1)
        #ax=DisplayArray(ax,hamiltonian)
        #plt.show()
        #exit()

    def reset_hartree(self):
        self.hartree=np.zeros([self.n_dof], dtype=complex_dtype)
        self.hartree_iterations=[]

    def reset_fock(self):
        self.fock=np.zeros([self.n_dof,self.n_dof], dtype=complex_dtype)
        self.fock_iterations=[]

    def reset_gorkov(self):
        self.gorkov=np.zeros([self.n_dof,self.n_dof], dtype=complex_dtype)
        self.gorkov_iterations=[]

    def set_hartree(self,onsite_amplitude,spins=None,orbitals=None,axes=None,location=None):
        self.hartree += np.diag(self._onsite(onsite_amplitude=onsite_amplitude, spin=spin, orbital=orbital, axes=axes, location=location))

    def set_fock(self, amplitude, orb_i=None, orb_f=None, hop_vector=None, add_time_reversal=True):
        self.fock += self._hopping(hopping_amplitude=amplitude, orb_i=orb_i, orb_f=orb_f, hop_vector=hop_vector, add_time_reversal=add_time_reversal)

    def set_gorkov(self, amplitude, orb_i=None, orb_f=None, hop_vector=None, add_time_reversal=True):
        self.gorkov += self._hopping(hopping_amplitude=amplitude, orb_i=orb_i, orb_f=orb_f, hop_vector=hop_vector, add_time_reversal=add_time_reversal)

    def set_external_hartree(self,onsite_amplitude,spins=None,orbitals=None,axes=None,location=None):
        self.external_hartree += np.diag(self._onsite(onsite_amplitude=onsite_amplitude, spin=spin, orbital=orbital, axes=axes, location=location))

    def set_external_fock(self, amplitude, orb_i=None, orb_f=None, hop_vector=None, add_time_reversal=True):
        self.external_fock += self._hopping(hopping_amplitude=amplitude, orb_i=orb_i, orb_f=orb_f, hop_vector=hop_vector, add_time_reversal=add_time_reversal)

    def set_external_gorkov(self, amplitude, orb_i=None, orb_f=None, hop_vector=None, add_time_reversal=True):
        self.external_gorkov += self._hopping(hopping_amplitude=amplitude, orb_i=orb_i, orb_f=orb_f, hop_vector=hop_vector, add_time_reversal=add_time_reversal)

    def set_hubbard_u(self, amplitude, orb_i=None, orb_f=None, hop_vector=None, add_time_reversal=True):
        if type(self._hubbard_u)==type(None):
            self._hubbard_u = np.zeros([self.n_dof, self.n_dof], dtype=complex_dtype) 
        self._hubbard_u += self._hopping(hopping_amplitude=amplitude, orb_i=orb_i, orb_f=orb_f, hop_vector=hop_vector, add_time_reversal=add_time_reversal)

    def _set_hubbard_indices(self):

        self._hubbard_indices = np.nonzero(self._hubbard_u)

        self.U_entries = self._hubbard_u[self._hubbard_indices]

        self.anomalous_indices=np.array(self._hubbard_indices)
        self.anomalous_indices+=self.n_dof
        self.anomalous_indices=tuple(self.anomalous_indices)

    def record_hartree(self, location, spin, orbital, _print=False):
        self.hartree_print = _print
        tmp=np.append(np.append(location,orbital),spin)
        index=np.ravel(coordinates_to_indices(tmp, self.extended_dimensions))
        self._hartree_indices.append(index)
        self._hartree_print_indices.append(_print)

    def record_fock(self, location_a, location_b, spin_a, spin_b, orbital_a=None, orbital_b=None, _print=False):
        self.fock_print = _print
        if type(orbital_a)==type(None) and type(orbital_b)==type(None):
            orbital_a=orbital_b=0
        tmp_a=np.append(np.append(location_a,orbital_a),spin_a)
        tmp_b=np.append(np.append(location_b,orbital_b),spin_b)
        index_a=coordinates_to_indices(tmp_a, self.extended_dimensions)
        index_b=coordinates_to_indices(tmp_b, self.extended_dimensions)
        temp=np.zeros([self.n_dof,self.n_dof])
        temp[index_a,index_b]=1
        temp = np.ravel(np.nonzero(temp))
        self._fock_indices.append(temp)
        self._fock_print_indices.append(_print)

    def record_gorkov(self, location_a, location_b, spin_a, spin_b, orbital_a=None, orbital_b=None, _print=False):
        self.gorkov_print = _print
        if type(orbital_a)==type(None) and type(orbital_b)==type(None):
            orbital_a=orbital_b=0
        tmp_a=np.append(np.append(location_a,orbital_a),spin_a)
        tmp_b=np.append(np.append(location_b,orbital_b),spin_b)
        index_a=coordinates_to_indices(tmp_a, self.extended_dimensions)
        index_b=coordinates_to_indices(tmp_b, self.extended_dimensions)
        temp=np.zeros([self.n_dof,self.n_dof])
        temp[index_a,index_b]=1
        temp = np.ravel(np.nonzero(temp))
        self._gorkov_indices.append(temp)
        self._gorkov_print_indices.append(_print)

    def set_friction(self,friction):
        self.friction = friction

    def set_max_iterations(self,max_iterations):
        self.max_iterations = max_iterations

    def set_absolute_convergence_factor(self,absolute_convergence_factor):
        self.absolute_convergence_factor = absolute_convergence_factor

    ########################################
    ######### Statistical Mechanics ########
    ########################################
        
    def _anomalous_density_matrix(self, v):
        return np.einsum('in,in->in',v[:self.n_dof],v[self.n_dof:])

    def _thermal_density_matrix(self, eigenvalues, eigenvectors: float):
        """Normal and anomalous thermal density matrices. Nonzero values
        are output as hartree, fock, gorkov fields.

        parameters
        ----------
        w : array
            eigenvalues
        dm : array
            wavefunction norm density matrix
        t : float
            temperature

        returns 
        ----------
        dm : array
            thermal density matrix
        """
        # Zero temperature thermal density matrix:
        # positive and negative eigenvalues are split
        # negative are selected
        # they're summed over
        w,v=eigenvalues,eigenvectors
        T=self.temperature
        tmp0=v[self._hubbard_indices[0]]
        tmp1=np.conj(v)[self._hubbard_indices[1]]
        tmp1a=v[self.anomalous_indices[1]]
        tmp0=tmp0[:,:self.n_dof]
        tmp1=tmp1[:,:self.n_dof]
        tmp1a=tmp1a[:,:self.n_dof]
        w=w[:self.n_dof]
        if T==0:
            trace_density = np.einsum('in,in->i',v[:self.n_dof,:self.n_dof],np.conj(v[:self.n_dof,:self.n_dof]),optimize=True)
            density = np.einsum('in,in->i',tmp0,tmp1,optimize=True)
            anomalous_density = -np.einsum('in,in->i',tmp0,tmp1a,optimize=True)

        # Nonzero temperature thermal density matrix:
        # the eigenvalues are weighted by the Fermi function
        else:
            # Fermi function:
            f = 1/(np.exp(w/T)+1)
            trace_density = np.einsum('in,in,n->i',v[:self.n_dof,:self.n_dof],np.conj(v[:self.n_dof,:self.n_dof]),f,optimize=True)
            density = np.einsum('in,in,n->i',tmp0,tmp1,f,optimize=True)
            anomalous_density = -np.einsum('in,in,n->i',tmp0,tmp1a,f,optimize=True)
        # print(trace_density)
        # print(density)
        # print(anomalous_density)
        # exit()
        return trace_density, density, anomalous_density

    def anomalous_fourier_transform(self,v):
        coeff=self.coeff
        coeff=np.block([[coeff,np.conj(coeff)],[np.conj(coeff),coeff]])
        return np.dot(coeff.T,v)

    def _set_fields(self, trace_density, density, anomalous_density):
        """Returns Hartree, Fock, Gorkov"""
        hartree = +np.einsum('ij,j->i',self._hubbard_u,trace_density,optimize=True)
        fock    = -np.multiply(self.U_entries,density)
        gorkov  = +np.multiply(self.U_entries,anomalous_density)

        return hartree, fock, gorkov
        # return np.real_if_close(hartree), np.real_if_close(fock), np.real_if_close(gorkov)

    def calculate_free_energy(self, w, trace_density, density, anomalous_density):

        self.Eg.append(np.real(np.sum(w[:self.n_dof] + np.diagonal(self._hamiltonian)[:self.n_dof])/self.n_dof))

        if self.temperature==0:
            self.f_mf=self.Eg
        else:
            self.f_mf=np.real(self.Eg-2*self.temperature*np.sum(np.log(1+np.exp(-w[self.n_dof:]/self.temperature)))/self.n_dof)
        
        h=np.sum(self.hartree*trace_density)
        h+=h
        f=np.sum(self.fock*density)
        f+=np.conj(f)
        g=np.sum(self.gorkov*anomalous_density)
        g+=np.conj(g)

        self.V_mf.append(np.real(-(h+f+g)/self.n_dof))

        h=np.einsum('ij,i,j->',self._hubbard_u,trace_density,trace_density,optimize=True)
        f=np.einsum('i,i,i->',self.U_entries,density,np.conj(density),optimize=True)
        g=np.einsum('i,i,i->',self.U_entries,anomalous_density,np.conj(anomalous_density),optimize=True)

        self.V.append(np.real(-(h-f+g)/self.n_dof))

        self.free_energy.append(self.f_mf[-1]+self.V[-1]-self.V_mf[-1])

        if self.print_V:
            print(self.V[-1])
        if self.print_V_mf:
            print(self.V_mf[-1])
        if self.print_Eg:
            print(self.Eg[-1])
        if self.print_free_energy:
            print(self.free_energy[-1])

        return self.free_energy[-1]

    def record_fields(self):

        hubbard_indices=np.transpose(self._hubbard_indices)

        temp_hartree=[]
        temp_fock=[]
        temp_gorkov=[]

        for index in self._hartree_indices:
            temp_hartree.append(complex(self.hartree[index]))
        
        for index in self._fock_indices:
            for i,hubbard_index in enumerate(hubbard_indices):
                if np.array_equal(index,hubbard_index):
                    temp_fock.append(self.fock[i])

        for index in self._gorkov_indices:
            for i,hubbard_index in enumerate(hubbard_indices):
                if np.array_equal(index,hubbard_index):
                    temp_gorkov.append(self.gorkov[i])

        self.hartree_iterations.append(temp_hartree)
        self.fock_iterations.append(temp_fock)
        self.gorkov_iterations.append(temp_gorkov)

        if np.sum(self._hartree_print_indices)>0:
            print('\nHartree field:')
            temp=np.array(self.hartree_iterations[-1])
            print(temp[self._hartree_print_indices])
        if np.sum(self._fock_print_indices)>0:
            print('\nFock field:')
            temp=np.array(self.fock_iterations[-1])
            print(temp[self._fock_print_indices])
        if np.sum(self._gorkov_print_indices)>0:
            print('\nGorkov field:')
            temp=np.array(self.gorkov_iterations[-1])
            print(temp[self._gorkov_print_indices])

    def __iter__(self):

        return self

    def __next__(self):
        
        ################## Record fields ###################
        
        self.record_fields()

        ############ Solve mean-field Hamiltonian #############

        self.set_mean_field_hamiltonian()
        w,v = self.solve()

        trace_density, density, anomalous_density = self._thermal_density_matrix(w,v)

        hartree, fock, gorkov = self._set_fields(trace_density, density, anomalous_density)

        ##################### Free energy ######################
        
        self.calculate_free_energy(w, trace_density, density, anomalous_density)

        ####################### Friction ########################

        hartree = (1-self.friction)*hartree + self.friction*self.hartree
        fock = (1-self.friction)*fock + self.friction*self.fock
        gorkov = (1-self.friction)*gorkov + self.friction*self.gorkov
        
        ###################### Update fields #####################

        self.hartree, self.fock, self.gorkov = hartree, fock, gorkov 

        return w,v

    def self_consistent_calculation(self):
        """If dos=True, the density of states are calculated once the self consistent
        loop has converged."""

        t = time.time()

        self._set_hubbard_indices()

        self.fock = self.fock[self._hubbard_indices]
        self.gorkov = self.gorkov[self._hubbard_indices]

        iteration = iter(self)
        
        for i in tqdm(range(self.max_iterations)):
            
            self.iterations+=1

            hartree_old = self.hartree
            fock_old = self.fock
            gorkov_old = self.gorkov

            w,v = next(iteration)
            
            eps = self.absolute_convergence_factor
            if (np.allclose(hartree_old,self.hartree,atol=eps) & np.allclose(fock_old,self.fock,atol=eps) & np.allclose(gorkov_old,self.gorkov,atol=eps)):
                self.converged=True
                break
            if i+1 == self.max_iterations:
                self.converged=False
                print('Did not converge within max_iterations!')
                break

        self.exec_time = time.time() - t

        self.hartree_iterations = np.transpose(self.hartree_iterations)
        self.fock_iterations = np.transpose(self.fock_iterations)
        self.gorkov_iterations = np.transpose(self.gorkov_iterations)

        return w,v

    def hartree(self):
        temp=np.reshape(self._hartree,self.extended_dimensions,'F')
        return np.real_if_close(temp)

    def fock(self):
        extended_dimensions=np.append(self.extended_dimensions,self.extended_dimensions)
        temp=np.zeros([self.n_dof,self.n_dof],dtype=complex_dtype)
        temp[self._hubbard_indices]=self._fock
        temp=np.reshape(temp,extended_dimensions,'F')
        return np.real_if_close(temp)

    def gorkov(self):
        extended_dimensions=np.append(self.extended_dimensions,self.extended_dimensions)
        temp=np.zeros([self.n_dof,self.n_dof],dtype=complex_dtype)
        temp[self._hubbard_indices]=self._gorkov
        temp=np.reshape(temp,extended_dimensions,'F')
        return np.real_if_close(temp)

class processed_data(bogoliubov_de_gennes):

    def __init__(self, model: "tight binding or bogoliubov_de_gennes", energy_interval: '1D array', resolution: complex, eigenvalues=None, eigenvectors=None, include_renormalisation_fields=True, include_k_space=True):
        """Calls LAPACK Hermitian matrix solver from scipy.linalg.eigh
        with overwrite to conserve memory

        Optional data attributes
        ----------

        density_of_states 
        ft_density_of_states 
        anomalous_density_of_states 
        ft_anomalous_density_of_states 
        hartree, fock, gorkov
        recorded_fields
        """
        
        self.dimensions = model.dimensions
        self.extended_dimensions = model.extended_dimensions
        self.basis = model.basis
        self.orbitals = model.orbitals
        self.hoppings = model.hoppings
        self.label_hoppings = model.label_hoppings
        self.impurities = model.impurities
        self.n_dof = model.n_dof
        self.n_dim = model.n_dim
        self.n_spins = model.n_spins
        self.n_orbitals = model.n_orbitals
        self.n_cells = model.n_cells
        self.coord = model.coord
        self.coord_cell = model.coord_cell
        self.coordinates_to_orbital_position = model.coordinates_to_orbital_position
        self.reciprocal_to_orbital_position = model.reciprocal_to_orbital_position
        self.centre = model.centre
        self.com = model.com

        self.hartree_iterations=model.hartree_iterations
        self.fock_iterations=model.fock_iterations
        self.gorkov_iterations=model.gorkov_iterations
        self.V=model.print_V
        self.V_mf=model.V_mf
        self.Eg=model.Eg
        self.free_energy=model.free_energy

        self.energy_interval=energy_interval
        self.resolution=resolution

        n_energy=len(energy_interval)

        if model._model=='bdg':
            
            self._hubbard_indices = model._hubbard_indices
            self.iterations = model.iterations

            if type(eigenvalues)!=type(None) and type(eigenvectors)!=type(None):

                anomalous_density_matrix = self._anomalous_density_matrix(eigenvectors) 
                ados = self._density_of_states(eigenvalues,anomalous_density_matrix)
                ados = np.real_if_close(ados)
                self.anomalous_density_of_states = np.reshape(ados, np.append(self.extended_dimensions,[n_energy]), 'F')

                if include_k_space:

                    ft_eigenvectors = model.fourier_transform_psi(eigenvectors)
                    ft_anomalous_density_matrix = self._anomalous_density_matrix(ft_eigenvectors) 
                    ft_ados = self._density_of_states(eigenvalues,ft_anomalous_density_matrix)
                    ft_ados = np.real_if_close(ft_ados)
                    self.ft_anomalous_density_of_states = np.reshape(ft_ados, np.append(self.extended_dimensions,[n_energy]), 'F')
                
            if include_renormalisation_fields:

                self._hartree_zip = np.real_if_close(model.hartree)
                self._fock_zip = np.real_if_close(model.fock)
                self._gorkov_zip = np.real_if_close(model.gorkov)
            
            # field recording:
            self.hartree_record = np.real_if_close(model.hartree_iterations)
            self.fock_record = np.real_if_close(model.fock_iterations)
            self.gorkov_record = np.real_if_close(model.gorkov_iterations)

        if type(eigenvalues)!=type(None) and type(eigenvectors)!=type(None):

            if include_k_space:

                ft_eigenvectors = model.fourier_transform_psi(eigenvectors)
                ft_density_matrix = self._density_matrix(ft_eigenvectors) 

                ft_dos = self._density_of_states(eigenvalues,ft_density_matrix)
                ft_dos = np.real_if_close(ft_dos)
                self.ft_density_of_states = np.reshape(ft_dos, np.append(self.extended_dimensions,[n_energy]), 'F')

            density_matrix = self._density_matrix(eigenvectors) 

            dos = self._density_of_states(eigenvalues,density_matrix)
            dos = np.real_if_close(dos)
            self.density_of_states = np.reshape(dos, np.append(self.extended_dimensions,[n_energy]), 'F')

    def hartree(self):
        temp = self._hartree_zip
        temp = np.reshape(temp, self.extended_dimensions, 'F')
        return temp

    def fock(self):
        temp=np.zeros([self.n_dof,self.n_dof],dtype=complex_dtype)
        temp[self._hubbard_indices]=self._fock_zip
        temp = np.reshape(temp, np.append(self.extended_dimensions,self.extended_dimensions), 'F')
        return temp

    def gorkov(self):
        temp=np.zeros([self.n_dof,self.n_dof],dtype=complex_dtype)
        temp[self._hubbard_indices]=self._gorkov_zip
        temp = np.reshape(temp, np.append(self.extended_dimensions,self.extended_dimensions), 'F')
        return temp

    def _greens_function(self, omega, eigenvalues, density_matrix):
        '''7-dimensional data set: [x, y, z, orbital, orbital, spin, spin]'''
        return np.einsum('e,ie->i', 1/(omega-eigenvalues), density_matrix,optimize=True)

    def _density_of_states(self, eigenvalues, density_matrix):
        '''8-dimensional data set: [x, y, z, orbital, orbital, spin, spin, omega]'''
        omegas = np.array(self.energy_interval, dtype=complex_dtype) + 1.0j*self.resolution
        green = np.array([self._greens_function(omega, eigenvalues, density_matrix) for omega in omegas])
        dos = -(1/np.pi)*np.imag(green)
        dos = np.moveaxis(dos,0,-1)
        return dos

    def local_density_of_states(self, energy, orbital, spin=None):
        temp = self.density_of_states
        if type(spin)==type(None):
            temp = np.sum(temp,-2)
        else:
            temp = temp[...,spin,:]
        index = FindNearestValueOfArray(self.energy_interval,energy)
        temp = temp[...,index]
        temp = temp[..., orbital]
        return temp

    def magnetism(self, energy, orbital):
        temp = self.density_of_states
        if self.n_spins==1:
            return temp[...,0,:]*0
        temp = temp[...,0,:]-temp[...,1,:]
        index = FindNearestValueOfArray(self.energy_interval,energy)
        temp = temp[...,index]
        temp = temp[..., orbital]
        return temp

    def spectrum(self, locations, orbital, spins=None):
        temp = self.density_of_states
        if type(spin)==type(None):
            temp = np.sum(temp,-2)
        else:
            temp = temp[...,spin,:]
        temp = temp[..., orbital,:]
        location=locations[0]
        temp = np.array([temp[(*location, ...)] for location in locations])
        return temp

    def ldos_cartesian(self, ldos, n_pts=40, gaussian_mean=0.4):

        def gaussian_convolution(x,x0,sigma):
            return np.exp(-np.sum([ (x[i]-x0[i])**2 for i in range(self.n_dim)],0)/(2*sigma**2))
        coord = self.coordinates_to_orbital_position
        xx = []
        n_pts=np.multiply([n_pts,n_pts],self.dimensions)
        for i in range(self.n_dim):
            _min = np.min(coord[...,i])-0.5
            _max = np.max(coord[...,i])+0.5
            xx.append(np.linspace(_min,_max,n_pts[i]))
        xx = np.meshgrid(*xx)

        coords = np.reshape(coord,[self.n_orbitals,self.n_cells,self.n_dim])
        ldos = np.reshape(ldos,[self.n_orbitals,self.n_cells])

        gaussian = [np.sum([gaussian_convolution(xx,coords[o,i],gaussian_mean)*ldos[o,i] for i in range(self.n_cells)],0) for o in range(self.n_orbitals)]
        
        gaussian=np.sum(gaussian,0)
        gaussian=np.flip(np.fft.fftshift(gaussian).T,axis=0)
        return gaussian

    def qpi_cartesian(self, ldos, n_pts=40, gaussian_mean=2, remove_central_bright_spot=True):
        
        ldos = np.fft.fft2(ldos, axes=(0,1), norm='ortho')

        if remove_central_bright_spot:
            ldos[0,0]=0

        def gaussian_convolution(x,x0,sigma):
            return np.exp(-np.sum([ (x[i]-x0[i])**2 for i in range(self.n_dim)],0)/(2*sigma**2))
        coord = self.reciprocal_to_orbital_position
        xx = []
        n_pts=np.multiply([n_pts,n_pts],self.dimensions)
        for i in range(self.n_dim):
            _min = np.min(coord[...,i])*1.1
            _max = np.max(coord[...,i])*1.1
            xx.append(np.linspace(_min,_max,n_pts[i]))
        xx = np.meshgrid(*xx)

        coords = np.reshape(coord,[self.n_orbitals,self.n_cells,self.n_dim])
        ldos = np.reshape(ldos,[self.n_orbitals,self.n_cells])

        gaussian = [np.sum([gaussian_convolution(xx,coords[o,i],gaussian_mean)*ldos[o,i] for i in range(self.n_cells)],0) for o in range(self.n_orbitals)]

        gaussian=np.sum(gaussian,0) #sum over orbitals

        gaussian = np.abs(gaussian)
        
        gaussian=np.sum(gaussian,-1)

        ldos=np.abs(ldos)

        return gaussian


class plot_data(processed_data):

    def __init__(self, fig, ax, data):
        
        self.fig = fig
        self.ax = ax
        
        self.data=data
        self.energy_interval=data.energy_interval
        self.resolution=data.resolution
        self.density_of_states=data.density_of_states
        self.n_spins = data.n_spins
        self.n_dim = data.n_dim
        self.hoppings = data.hoppings
        self.label_hoppings = data.label_hoppings
        self.impurities = data.impurities
        self.n_orbitals = data.n_orbitals
        self.orbitals = data.orbitals
        self.n_cells = data.n_cells
        self.dimensions = data.dimensions
        self.coord = data.coord
        self.coord_cell = data.coord_cell
        self.centre = data.centre
        self.coordinates_to_orbital_position = data.coordinates_to_orbital_position
        self.basis = data.basis
        self.com = data.com

        self.hartree_iterations=data.hartree_iterations
        self.fock_iterations=data.fock_iterations
        self.gorkov_iterations=data.gorkov_iterations
        self.V=data.V
        self.V_mf=data.V_mf
        self.Eg=data.Eg
        self.free_energy=data.free_energy

    def _imshow(self, ldos):

        self.interpolation = 'none'

        x=np.fft.fftshift(ldos)

        self.im = self.ax.imshow(
                x[::-1].T, extent=self._extent,
                interpolation=self.interpolation,
                vmin=self.vmin, vmax=self.vmax,
                cmap=self.cmap,
                origin='lower')

    def differential_current_map(self, energy, layer=(ALL,ALL,0), orbital=0, spin=None, cartesian=False, n_pts=10, gaussian_mean=0.4):

        x,y=self.dimensions[:2]
        self._extent = [-x/2,x/2,-y/2,y/2]

        ldos = self.data.local_density_of_states(energy, orbital, spin)

        ldos = ldos[layer]
        if cartesian:
            ldos=self.data.ldos_cartesian(ldos, n_pts, gaussian_mean)
            self._extent = [-x/2-0.5/n_pts,x/2+0.5/n_pts,-y/2-0.5/n_pts,y/2+0.5/n_pts]

        eV = energy
        epsilon = self.resolution

        self.text= (f'DOS map'
        '\n'
        f'$\omega={eV:.2f}$'
        '\n'
        f'$\epsilon={epsilon}$'
        )    

        self.vmin, self.vmax = np.amin(ldos), np.amax(ldos)

        if layer==(ALL,ALL,0):
            self.ax.set(xlabel='$x/a_x$', ylabel='$y/a_y$')

        if layer==(ALL,0,ALL):
            self.ax.set(xlabel='$x/a_x$', ylabel='$z/a_z$')

        if layer==(0,ALL,ALL):
            self.ax.set(xlabel='$y/a_y$', ylabel='$z/a_z$')

        self.cmap = ListedColormap(cm.afmhot(np.linspace(0, 0.75, 256)))

        self._imshow(ldos)

        self.fig.subplots_adjust(hspace=0.1, wspace=0.1)

        return self.fig,self.ax

    def probe_magnetic_z(self, energy, layer=(ALL,ALL,0), orbital=0):

        ldos = self.data.local_density_of_states(energy, orbital)

        ldos = ldos[...,0] - ldos[...,1]
    
        ldos = ldos[layer]

        eV = energy
        epsilon = self.resolution

        self.text= (f'\map'
        '\n'
        f'$\omega={eV:.2f}$'
        '\n'
        f'$\epsilon={epsilon}$'
        )    

        self.vmin, self.vmax = np.amin(ldos), np.amax(ldos)

        x,y=np.shape(ldos)[:2]

        if layer==(ALL,ALL,0):
            self.ax.set(xlabel='$x/a_x$', ylabel='$y/a_y$')

        if layer==(ALL,0,ALL):
            self.ax.set(xlabel='$x/a_x$', ylabel='$z/a_z$')

        if layer==(0,ALL,ALL):
            self.ax.set(xlabel='$y/a_y$', ylabel='$z/a_z$')

        self._extent = [-x/2,x/2,-y/2,y/2]
        
        self.cmap = ListedColormap(cm.PiYG(np.linspace(0, 1, 256)))

        self._imshow(ldos)

        self.fig.subplots_adjust(hspace=0.1, wspace=0.1)

        return self.fig,self.ax

    def quasiparticle_interference(self, energy, layer=(ALL,ALL,0), orbital=0, spins=None, remove_central_bright_spot=True):
        
        ldos = self.ldos=greens_function.local_density_of_states(energy, orbital, spin)
        ldos = self.ldos[layer]

        self.vmin, self.vmax = np.amin(ldos), np.amax(ldos)

        f = np.fft.fft2(ldos, axes=(0,1), norm='ortho')
        abs_f = np.abs(f)

        if layer==(ALL,ALL,0):
            self.ax.set(xlabel='$k_x a_x$', ylabel='$k_y a_y$')

        if layer==(ALL,0,ALL):
            self.ax.set(xlabel='$k_x a_x$', ylabel='$k_z a_z$')

        if layer==(0,ALL,ALL):
            self.ax.set(xlabel='$k_y a_y$', ylabel='$k_z a_z$')
        
        self._extent = [-1,1,-1,1]

        x_label_list = ['$-\pi$', '0', '$\pi$']
        y_label_list = ['$-\pi$', '0', '$\pi$']

        self.ax.set_xticks([-1,0,1])
        self.ax.set_yticks([-1,0,1])

        self.ax.set_xticklabels(x_label_list)
        self.ax.set_yticklabels(x_label_list)


        # max/min without central bright spot and lines:
        if remove_central_bright_spot:
            vmax=[]
            vmin=[]
            temp=np.copy(abs_f)
            temp[0,0]=0
            vmax.append(np.amax(temp))
            vmin.append(np.amin(temp))
            self.vmax = np.max(vmax)
            self.vmin = np.max(vmin)
        
        self.cmap = LinearSegmentedColormap.from_list("", ["black","orange","white"])

        self._imshow(abs_f)

        self.fig.subplots_adjust(hspace=0.1, wspace=0.1)

        return self.fig,self.ax
        
    def reciprocal_space_surface(self, energy, layer=(ALL,ALL,0), orbital=0, spins=None):

        ldos = self.ldos=greens_function.local_density_of_states(energy, orbital, spin)

        ldos = self.ldos[layer]

        self.vmin, self.vmax = np.amin(ldos), np.amax(ldos)

        x,y=np.shape(ldos)[:2]

        if layer==(ALL,ALL,0):
            self.ax.set(xlabel='$k_x a_x$', ylabel='$k_y a_y$')

        if layer==(ALL,0,ALL):
            self.ax.set(xlabel='$k_x a_x$', ylabel='$k_z a_z$')

        if layer==(0,ALL,ALL):
            self.ax.set(xlabel='$k_y a_y$', ylabel='$k_z a_z$')

        self._extent = [-1,1,-1,1]

        x_label_list = ['$-\pi$', '0', '$\pi$']
        y_label_list = ['$-\pi$', '0', '$\pi$']

        self.ax.set_xticks([-1,0,1])
        self.ax.set_yticks([-1,0,1])

        self.ax.set_xticklabels(x_label_list)
        self.ax.set_yticklabels(x_label_list)
        
        self.cmap = ListedColormap(cm.afmhot(np.linspace(0, 0.75, 256)))

        self._imshow(ldos)

        self.fig.subplots_adjust(hspace=0.1, wspace=0.1)

        return self.fig,self.ax

    def band_structure(self, dimension, oribital=0, spins=None):

        self.cmap = ListedColormap(cm.afmhot(np.linspace(0, 0.75, 256)))
        self.interpolation = 'none'

        self.xlabel=['$x/a_x$','$y/a_y$','$z/a_z$'][dimension]

        self.ylabel='$\omega$'
        
        dos = self.data.density_of_states
        self.extended_dimensions=np.arange(len(np.shape(dos))-1)
        self.extended_dimensions=tuple(np.delete(self.extended_dimensions,dimension))
        dos=np.sum(dos,axis=self.extended_dimensions)

        x=np.shape(dos)[0]
        y=np.real(self.energy_interval)

        self._extent = [-x/2,x/2,y[0],y[-1]]

        x=dos
        x=np.fft.fftshift(x,0)
        self.im = self.ax.imshow(
                x.T, extent=self._extent, origin='lower',
                interpolation=self.interpolation,
                #vmin=self.vmin, vmax=self.vmax,
                cmap=self.cmap)

        self.ax.set(xlabel=self.xlabel, ylabel=self.ylabel)

        return self.fig, self.ax

    def set_text_box(self):
        text_box = AnchoredText(self.text, loc=2, pad=0.3, borderpad=0)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        self.ax.add_artist(text_box)

    def set_cbar(self): 
        cax = inset_axes(self.ax,width = '5%',height = '80%',loc = 5) 
        self.cbar = plt.colorbar(self.im, ax = self.ax,cax = cax ,format = '%.2f' ,ticks = [self.vmin,self.vmax])       
        cax.yaxis.set_ticks_position('left')

    def set_label(self,label):
        text_box = AnchoredText(label, loc=3, pad=0.3, borderpad=0, frameon=False)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        self.ax.add_artist(text_box)

    def plot_3D(self, energy, orbital=0, spins=None):

        import plotly.graph_objects as go

        ldos = self.ldos=greens_function.local_density_of_states(energy, orbital, spin)

        x=int(self.dimensions[0]/2)
        y=int(self.dimensions[1]/2)
        z=int(self.dimensions[2]/2)
        X, Y, Z = np.mgrid[-x:x+1, -y:y+1, -z:z+1]
        values = np.fft.fftshift(self.ldos)
        #values=self.ldos
        minv=np.min(values)
        maxv=np.max(values)

        fig = go.Figure(data=go.Volume(
            x=X.flatten(),
            y=Y.flatten(),
            z=Z.flatten(),
            value=values.flatten(),
            isomin=maxv/2,
            isomax=maxv,
            opacity=1, # needs to be small to see through all surfaces
            opacityscale='extremes',
            surface_count=10, # needs to be a large number for good volume rendering
            caps= dict(x_show=False, y_show=False, z_show=False), # no caps
            ))

        fig.update_layout(
            scene = dict(
                xaxis = dict(nticks=3, tickvals=[-x,0,x],),
                yaxis = dict(nticks=3, tickvals=[-y,0,y],),
                zaxis = dict(nticks=3, tickvals=[-z,0,z],),))
            #width=700,
            #margin=dict(r=20, l=10, b=10, t=10))

        return fig

    def fields(self, xaxis, xlabel, field, label, twin_field=None, twin_label=None, second_twin_field=None, second_twin_label=None):
    
        color='tab:red'
        
        if len(field)!=len(label):
            raise ValueError("length of field not equal to length of label")

            if type(xaxis)==type(None):
                xaxis = np.arange(len(field))

            field=field
            label=label
            self.ax.plot(xaxis, field, color=color,marker='.',markersize=4,label=label)
            self.ax.set_ylabel(label, color=color)
            self.ax.tick_params(axis='y', labelcolor=color)
        
        if type(twin_field)!=type(None):

            self.twin_ax = self.ax.twinx()
            
            color = 'tab:blue'

            if len(twin_field)!=len(twin_label):
                raise ValueError("length of twin_field not equal to length of twin_label")
            
            field=twin_field
            label=twin_label

            self.twin_ax.plot(xaxis, field, color=color,marker='x',markersize=4,label=label)  
            self.twin_ax.set_ylabel(label, color=color)
            self.twin_ax.tick_params(axis='y', labelcolor=color)

        if type(second_twin_field)!=type(None):

            self.twin_ax2 = self.ax.twinx()
            
            color = 'tab:green'

            if len(second_twin_field)!=len(second_twin_label):
                raise ValueError("length of second_twin_field not equal to length of second_twin_label")

            field=second_twin_field
            label=second_twin_label

            self.twin_ax2.plot(xaxis, field, color=color,marker='x',markersize=4,label=label)  
            self.twin_ax2.set_ylabel(label, color=color)
            self.twin_ax2.tick_params(axis='y', labelcolor=color)

            self.twin_ax2.spines.right.set_position(("axes", 1.2))


        self.ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        self.ax.set_xlabel(xlabel)

        plt.tight_layout()

        return self.fig, self.ax 

    def spectrum(self, locations, orbital, spins=None):

        spectrum = self.data.spectrum(locations, orbital, spin)

        n = len(spectrum)
        colors = [cm.nipy_spectral(i) for i in np.linspace(0, 1, n)]

        for i in range(n):


            if not spin:
                s=0 
                x,y=self.energy_interval,spectrum[i,s]
                self.ax.plot(x,y, linestyle='solid', color=colors[i],label=(f'$\mathbf{{r}}={locations[i]}$, ' + r'$\uparrow$'))
                self.ax.fill_between(x,y,0,color=colors[i], alpha=0.4)
                s=1
                x,y=self.energy_interval,spectrum[i,s]
                self.ax.plot(x,y, linestyle='dotted', color=colors[i],label=(f'$\mathbf{{r}}={locations[i]}$, ' + r'$\downarrow$'))
                self.ax.fill_between(x,y,0,color=colors[i], alpha=0.4)
            else:
                x,y=self.energy_interval,spectrum[i]
                self.ax.plot(x,y, linestyle='solid', color=colors[i],label=(f'$\mathbf{{r}}={locations[i]}$'))
                self.ax.fill_between(x,y,0,color=colors[i], alpha=0.4)


        self.ax.set(xlabel=r'$\omega$')  
        self.ax.set(ylabel=r'Density of states')

        text_DOS = ('Spectrum'
        '\n'
        f'$\epsilon={self.resolution}$'
        )    

        self.fig.text(.22, .81, text_DOS,
                 {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
                               ec="none", pad=0.2)}, ha='center', va='center')
                               
        legend = self.fig.legend(loc="upper left", 
        fancybox=True, shadow=True, #prop=fontP,
        bbox_to_anchor=(0.14,0.58,1,0.2))

        return self.fig,self.ax

    def lattice(self, energy, model, layer, spins, orbitals, annotate_orbs=True, show_cell_borders=True, show_basis_vectors=True, show_hopping=True, show_impurities=True, show_only_centre=True):

        ldos = self.data.local_density_of_states(energy, orbitals)
        ldos = ldos[layer]

        # Colors and labels:
        self.color_orb = plt.cm.rainbow(np.linspace(0,1,self.n_orbitals))
        self.label_basis = list(string.ascii_uppercase)
        self.label_orb = list(string .ascii_lowercase)

        # Show hopping:
        if show_hopping:
            for i, link in enumerate(self.hoppings):
                self._plt_hopping(*link, self.label_hoppings[i], show_only_centre)

        # Find orbitals in layer to be plotted:
        for o, orb in enumerate(self.orbitals):
            tmp=layer+(o,slice(2))
            print(np.shape(self.coordinates_to_orbital_position))
            temp = self.coordinates_to_orbital_position[tmp]
            print(np.shape(temp))
            pos = temp.reshape([self.dimensions[0]*self.dimensions[1], 2])
            size = ldos[...,o].reshape([self.dimensions[0]*self.dimensions[1]])
            xx, yy = zip(*pos)
            s=size*40/(np.mean(size))
            self.ax.scatter(xx, yy, s=s, color=self.color_orb[o])

        # Impurities
        if show_impurities:
            for impurities in self.impurities:
                [impurity_loc, impurity_spin, impurity_orb] = impurities
                ll=[]
                for loc in impurity_loc:
                    if len(loc)==2:
                        loc.append(0)
                    ll.append(loc)
                ll = np.array(ll)
                oo=[]
                for orb in impurity_orb:
                    oo.append(orb)
                n_imp = np.shape(ll)[0]
                for o in oo:
                    if ll!=[]:
                        indices = tuple([list(ll[:,i]) for i in range(self.n_dim)]+[o,slice(None,2,1)])
                        pos = self.coordinates_to_orbital_position[indices]
                        xx, yy = zip(*pos)
                        s=10/np.mean(self.dimensions[0:1])
                        self.ax.scatter(xx, yy, s=6*s, color='k', marker='o')
                        self.ax.scatter(xx, yy, s=s, color=self.color_orb[o], marker='*')
        
        # Cell borders:
        if show_cell_borders==True:
            width = self.basis[0][0]
            # width = 2*basis[0][0] # 2 for hexagonal self
            height = self.basis[1][1]

            tmp=layer+(0,slice(2))
            temp = self.coordinates_to_orbital_position[tmp]
            for a_x, a_y in zip(*(temp.reshape([self.dimensions[0]*self.dimensions[1],2])+self.com[:2]).T):
                self.ax.add_patch(Rectangle(
                    xy=(a_x-width/2, a_y-height/2) ,width=width, height=height,
                    linewidth=1, color='blue', fill=False))
        
        # Annotate orbitals:
        for o, orb in enumerate(self.orbitals):
            if annotate_orbs==True:
                self.ax.annotate(self.label_orb[o], self.orbitals[o][:2]+np.array([0.02,0.02]), horizontalalignment='right', verticalalignment='bottom', color=self.color_orb[o])

        # Draw basis vectors:
        if show_basis_vectors==True:
            for i in range(2): # Two basis vectors are plotted in the plane
                self.ax.annotate(text='', xytext=self.basis[i][:2],xy=self.orbitals[0][:2],  arrowprops=dict(arrowstyle='<-', lw=2))
                self.ax.annotate(text=self.label_basis[i], xytext=0.5*self.basis[i][:2]+[0.05,0.05],xy=self.orbitals[0][:2])

        return self.fig, self.ax

    def _plt_hopping(self, orb_i: int, orb_f: int, cell_hop: tuple, label='', PBC=True, show_only_centre=True):
        cell_hop = np.array(cell_hop)
        cells = self.n_cells
        if (cell_hop==np.zeros([self.n_dim])).all() and PBC==True:
            if show_only_centre:
                index=[0]
            else:
                index=np.arange(cells)
            for i in index:
                index_i = self.coord[i][:self.n_dim]
                index_f = np.mod(index_i + cell_hop, self.dimensions)
                index_i = list(index_i)+[orb_i]
                orb_0 = self.coordinates_to_orbital_position[tuple(index_i)] 
                index_f = list(index_f)+[orb_f]
                orb_1 = self.coordinates_to_orbital_position[tuple(index_f)]
                orb_0=orb_0[:2]
                orb_1=orb_1[:2]
                self.ax.annotate(text='', xytext=orb_1, xy=orb_0,  arrowprops={'arrowstyle': '<->', 'ls': 'dashed'})
                self.ax.annotate(text=label, xytext=0.5*(orb_0+orb_1)+(0,0.01), xy=orb_0)
        # else:
        for dimension in range(self.n_dim):
            temp=np.zeros(self.n_dim)
            temp[dimension]=1
            if (temp==cell_hop).all():
                if show_only_centre:
                    index=[0]
                else:
                    index=np.arange(cells)
                for i in index: 
                    index_i = self.coord[i][:self.n_dim]
                    index_f = np.mod(index_i + cell_hop, self.dimensions)
                    coord_i = self.coord_cell[i]
                    edge = bool(coord_i[dimension]>=self.centre[dimension])
                    if edge and not PBC:
                        pass
                    else:
                        index_i = list(index_i)+[orb_i,ALL]
                        orb_0 = self.coordinates_to_orbital_position[tuple(index_i)] 
                        index_f = list(index_f)+[orb_f,ALL]
                        orb_1 = self.coordinates_to_orbital_position[tuple(index_f)]
                        orb_0=orb_0[:2]
                        orb_1=orb_1[:2]
                        self.ax.annotate(text='', xytext=orb_1, xy=orb_0,  arrowprops={'arrowstyle': '<->', 'ls': 'dashed'})
                        self.ax.annotate(text=label, xytext=0.5*(orb_0+orb_1)+(0,0.01), xy=orb_0)

