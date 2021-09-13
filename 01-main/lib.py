#!/usr/bin/env python

# PythU python self-consistent Hubbard U module
# September 12th, 2021
__version__='1.0.0'
# 
# To add in next version: 
# -function of momentum for U, mean-fields and tight binding
# -function of space for placing impurities
# -plotting classes

# Simulation of non-unitary, spin-triplet unconventional 
# superconductivity and magnetism for local density of states,
# quasiparticle interference and topological features.
#
# Copyright (C) 2021, Henry Joseph Sheehy
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
########################## Preconfig ##############################
###################################################################
Pauli_z = np.array([[1,0],[0,-1]])
Pauli_y = np.array([[0,-1.0j],[1.0j,0]])
Pauli_x = np.array([[0,1],[1,0]])
Pauli_plus = (Pauli_x+1.0j*Pauli_y)*(1./2.)
Pauli_minus = (Pauli_x-1.0j*Pauli_y)*(1./2.)
Pauli_vec = np.dstack([Pauli_x, Pauli_y, Pauli_z])
Pauli_vec = np.moveaxis(Pauli_vec,-1,0)
###################################################################
############################ General ##############################
###################################################################
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

def kron(a : 'array', b : 'array'):
    return np.kron(a,b)
    # return np.multiply.outer(a,b).reshape(np.multiply(np.shape(a),np.shape(b)))

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
##################### Statistical Mechanics #######################
###################################################################
def Green_function(omega, eigenvalues, density_matrix):
    '''7-dimensional data set: [x, y, z, spin, spin, orbital, orbital]'''
    return np.einsum('e,ie->i', 1/(omega-eigenvalues), density_matrix,optimize=True)

def DOS(energy_interval, resolution, eigenvalues, density_matrix):
    '''8-dimensional data set: [x, y, z, spin, spin, orbital, orbital, omega]'''
    omegas = np.array(energy_interval, dtype=complex_dtype) + 1.0j*resolution
    green = np.array([Green_function(omega, eigenvalues, density_matrix) for omega in omegas])
    dos = -(1/np.pi)*np.imag(green)
    dos = np.moveaxis(dos,0,-1)
    return dos

def local_density_of_states(density_of_states, index, trace_over=True):
    """[Spin-orbital resolved] local density of states
    
    parameters
    ----------
    omega : float
        finds nearest value of omega in omegas
    traced : bool
        if True traces over the spin-orbital componenets 
    layer : int
        z component

    returns 
    ----------
    ldos : array
        local density of states
    """
    ldos = density_of_states[...,index]
    if trace_over==True:
        tmp = np.shape(ldos)
        n_spins, n_orbs = tmp[-2:]
        ldos = np.sum(ldos, axis=(-1,-2))/(n_spins*n_orbs)
    return ldos

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
############################# Lattice #############################
###################################################################
#class tight_binding():
    """A class for the lattice containing the following objects:
    n_orbs
    n_dim
    n_cells
    n_atoms
    n_dof
    index: initialised Zn_Z
    coord: initialised Z_Zn
    coord_cell: array of coordinates in basis
    coord_orb: array of indices of orbitals
    pos_cell: array of real space coordinates of cells
    pos_orb: array of real space coordinates of orbitals
    com_orb: centre of mass of the orbitals in real space
    color_orb: list of colors 
    label_orb: list of strings
    label_basis: list of strings"""
#    def __init__(self, dimensions, n_spins, basis, orbitals, hoppings, impurities,pbc=None):

        #positions = np.empty(np.append(self.extended_dimensions,[3]))
        #self.pos_orb = np.copy(positions)
        #for index in range(self.n_dof):
        #    coord = tuple(self.coord[index])
            # Add zero vector to 2D basis:
        #    pos = np.dot(self.coord_cell[index], self.basis)
        #    positions[coord] = pos
        #    self.pos_orb[coord] = pos + self.orbitals[coord[4]]
        # Origin is the lattice centre of mass:
        #self.com = np.mean(orbitals,axis=0)
###################################################################
######################### Hamiltonian #############################
###################################################################
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

    def __init__(self, dimensions, n_spins, basis=None, orbitals=None, pbc=None): #, k_dim=None):
    
        self.dimensions=dimensions
        self.n_spins=n_spins
        self.basis=np.array(basis,dtype=float)
        self.orbitals=orbitals
        self.n_orbs = len(orbitals)
        self.n_dim = len(dimensions)
        self.temperature = 0 

        self.reciprocal_basis = la.inv(self.basis)

        #self.k_dim = k_dim

        if type(pbc)==type(None):
            self.pbc=np.ones([self.n_dim],dtype=bool)
        else:
            self.pbc=pbc

        print(self.pbc)

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
        self.extended_dimensions  = np.append(self.dimensions, [self.n_spins,self.n_orbs])
        self.n_cells = np.prod(self.dimensions)
        self.n_atoms = self.n_cells * self.n_orbs
        self.n_dof = self.n_atoms * n_spins

        self.centre = np.array(np.array(self.dimensions)/2, dtype=int)  # centred-coordinates
        self.coord=np.reshape(np.indices(self.extended_dimensions),[len(self.extended_dimensions),self.n_dof],'F').T

        self.coord[:,:self.n_dim] = wrap(self.coord[:,:self.n_dim], self.dimensions)

        self.coord_cell = self.coord[:,:self.n_dim]
        
        self.position_cell = np.dot(self.coord_cell[:self.n_cells],basis)

        self.position = np.array([[self.position_cell[i] + self.orbitals[m] for m in range(self.n_orbs)] for i in range(self.n_cells)])

        self.coordinates_to_cell_position = np.reshape(self.position_cell, [*self.dimensions,self.n_dim],'F')
        self.coordinates_to_orbital_position = np.reshape(self.position, [*self.dimensions,self.n_orbs,self.n_dim],'F')
        
        ####################################################

        self._hamiltonian = None

        self.coeff = self._fourier_coeffs()

    def set_temperature(self,temperature):
        self.temperature=temperature

    #################################### 
    ########### Hamiltonian ############
    ####################################

    def _onsite(self, onsite_amplitude, spin=None, orbital=None, axes=None, location=None):

        if np.isscalar(onsite_amplitude):
            if type(orbital)==type(None):
                orbit_tensor = np.eye(self.n_orbs)
            else:
                orbit_tensor = np.zeros([self.n_orbs,self.n_orbs])
                orbit_tensor[orbital]=1
            if type(spin)==type(None):
                spin_tensor = np.eye(self.n_spins)
            else:
                spin_tensor = np.zeros([self.n_spins,self.n_spins])
                spin_tensor[spin]=1
            temp = onsite_amplitude*np.kron(spin_tensor,orbit_tensor)

        elif np.shape(onsite_amplitude)[0]==self.n_spins:
            spin_tensor = hopping_amplitude
            if self.n_orbs>1:
                if type(orbital)==type(None):
                    ValueError('orbtal or spin_orbit_tensor not given!')
                else:
                    orbit_tensor=np.zeros([self.n_orbs, self.n_orbs])
                    orbit_tensor[orbital]=1
                    temp = kron(spin_tensor,orbit_tensor)
            else:
                temp = spin_tensor
        elif np.shape(hopping_amplitude)[0]==self.n_spins*self.n_orbs:
            temp = hopping_amplitude

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

        return np.kron(temp, sites)

    def set_onsite(self, onsite_amplitude, spin=None, orbital=None, axes=None, location=None):
        if type(self._hamiltonian)==type(None):
            self._hamiltonian = np.zeros([self.n_dof, self.n_dof], dtype=complex_dtype) 
        self._hamiltonian += self._onsite(onsite_amplitude, spin, orbital, axes, location)
        
    
    def set_hopping(self, hopping_amplitude, orb_i=None, orb_f=None, hop_vector=None, add_time_reversal=True):
        if type(self._hamiltonian)==type(None):
            self._hamiltonian = np.zeros([self.n_dof, self.n_dof], dtype=complex_dtype) 
        self._hamiltonian += self._hopping(hopping_amplitude=hopping_amplitude, orb_i=orb_i, orb_f=orb_f, hop_vector=hop_vector, add_time_reversal=add_time_reversal)
        
    def _hopping(self, hopping_amplitude, orb_i=None, orb_f=None, hop_vector=None, add_time_reversal=True):
            
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
            spin_tensor = hopping_amplitude
            if self.n_orbs>1:
                if type(orb_i)==type(None) or type(orb_f)==type(None):
                    ValueError('orb_i, orb_f or spin_orbit_tensor not given!')
                else:
                    orbit_tensor=np.zeros([self.n_orbs, self.n_orbs])
                    orbit_tensor[orb_i,orb_f]=1
                    spin_orbit_tensor=kron(spin_tensor,orbit_tensor)
            else:
                spin_orbit_tensor = spin_tensor
        elif np.shape(hopping_amplitude)[0]==self.n_spins*self.n_orbs:
            spin_orbit_tensor = hopping_amplitude

        temp=kron(spin_orbit_tensor,temp)
        
        onsite = bool(orb_i==orb_f and np.array_equal(hop_vector, np.zeros([self.n_dim])))
        if add_time_reversal and not onsite:
            temp+=dagger(temp)

        return temp

    def set_impurities(self, impurity_amplitude, locations, spin=None, orbital=None):
        for location in locations:
            self.set_onsite(onsite_amplitude=impurity_amplitude, spin=spin, orbital=orbital, location=location)

    def set_magnetic_impurities(self, M, locations, spin_matrix=None, orbital=None):
        spin_matrix=np.einsum('i,ijk->jk',M,Pauli_vec)
        self.set_onsite_spin(spin_matrix, orbital, locations)

    def unravel_hamiltonian(self):
        self._hamiltonian = np.reshape(self._hamiltonian, np.append(self.extended_dimensions,self.extended_dimensions), 'F')

    def ravel_hamiltonian(self):
        self._hamiltonian = np.reshape(self._hamiltonian, [self.n_dof,self.n_dof], 'F')

    def _fourier_coeffs(self):
        x=self.coord_cell
        k=np.copy(x)
        k=k*2*np.pi/(self.dimensions)-np.pi
        c=np.dot(x,k.T)
        c=np.exp(-1.0j*c)
        return c

    def fourier_transform_psi(self,v):
        coeff=self.coeff
        return np.dot(np.conj(coeff), v)

    def inv_fourier_transform_psi(self,v):
        coeff=self.coeff
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
        # v=self.fourier_transform_psi(v)
        # v=self.inv_fourier_transform_psi(v)
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
        
        self.friction = 0
        self.max_iterations = 100
        self.absolute_convergence_factor = 0.001
        
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

    def set_mean_field_hamiltonian(self):

        # initialise tight binding hamiltonian if doesn't exist:
        try: 
            self._tb_ham
        except:
            self._tb_ham = self._hamiltonian

        n_dof = self.n_dof

        try:
            self.hubbard_indices
        except:
            self._hubbard_u = np.diag(self.hartree)+self.fock+self.gorkov
            self._set_hubbard_indices()
            self.fock=self.fock[self.hubbard_indices]
            self.gorkov=self.gorkov[self.hubbard_indices]

        hartree=self.hartree
        fock=self.fock
        gorkov=self.gorkov

        hamiltonian = np.zeros([2*n_dof,2*n_dof],dtype=complex_dtype)
        hamiltonian[:n_dof,:n_dof]=self._tb_ham+np.diag(hartree)
        hamiltonian[n_dof:,n_dof:]=-np.conj(self._tb_ham)-np.conj(np.diag(hartree))
        hamiltonian[self.hubbard_indices[0],self.anomalous_indices[1]]=fock
        hamiltonian[self.anomalous_indices[0],self.hubbard_indices[1]]=-np.conj(fock)
        hamiltonian[self.hubbard_indices[0],self.anomalous_indices[1]]=-np.conj(gorkov)
        hamiltonian[self.anomalous_indices[0],self.hubbard_indices[1]]=gorkov
        self._hamiltonian = hamiltonian

    def set_hartree(self,onsite_amplitude,spin=None,orbital=None,axes=None,location=None):
        self.hartree += np.diag(self._onsite(onsite_amplitude=onsite_amplitude, spin=spin, orbital=orbital, axes=axes, location=location))

    def set_fock(self, amplitude, orb_i=None, orb_f=None, hop_vector=None, add_time_reversal=True):
        self.fock += self._hopping(hopping_amplitude=amplitude, orb_i=orb_i, orb_f=orb_f, hop_vector=hop_vector, add_time_reversal=add_time_reversal)

    def set_gorkov(self, amplitude, orb_i=None, orb_f=None, hop_vector=None, add_time_reversal=True):
        self.gorkov += self._hopping(hopping_amplitude=amplitude, orb_i=orb_i, orb_f=orb_f, hop_vector=hop_vector, add_time_reversal=add_time_reversal)

    def set_external_hartree(self,onsite_amplitude,spin=None,orbital=None,axes=None,location=None):
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

        self.hubbard_indices = np.nonzero(self._hubbard_u)

        self.U_entries = self._hubbard_u[self.hubbard_indices]

        self.anomalous_indices=np.array(self.hubbard_indices)
        self.anomalous_indices+=self.n_dof
        self.anomalous_indices=tuple(self.anomalous_indices)

    def record_hartree(self, location, spin, orbit, _print=False):
        self.hartree_print = _print
        tmp=np.append(np.append(location,spin),orbit)
        index=np.ravel(coordinates_to_indices(tmp, self.extended_dimensions))
        self._hartree_indices.append(index)
        self._hartree_print_indices.append(_print)

    def record_fock(self, location_a, location_b, spin_a, spin_b, orbital_a=None, orbital_b=None, _print=False):
        self.fock_print = _print
        if type(orbital_a)==type(None) and type(orbital_b)==type(None):
            orbital_a=orbital_b=0
        tmp_a=np.append(np.append(location_a,spin_a),orbital_a)
        tmp_b=np.append(np.append(location_b,spin_b),orbital_b)
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
        tmp_a=np.append(np.append(location_a,spin_a),orbital_a)
        tmp_b=np.append(np.append(location_b,spin_b),orbital_b)
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
        tmp0=v[self.hubbard_indices[0]]
        tmp1=np.conj(v)[self.hubbard_indices[1]]
        tmp1a=v[self.anomalous_indices[1]]
        if T==0:
            tmp0=tmp0[:,:self.n_dof]
            tmp1=tmp1[:,:self.n_dof]
            tmp1a=tmp1a[:,:self.n_dof]
            hartree_entries = np.einsum('in,in->i',v[:self.n_dof,:self.n_dof],np.conj(v[:self.n_dof,:self.n_dof]),optimize=True)
            fock_entries = -np.einsum('in,in->i',tmp0,tmp1,optimize=True)
            gorkov_entries = np.einsum('in,in->i',tmp0,tmp1a,optimize=True)

        # Nonzero temperature thermal density matrix:
        # the eigenvalues are weighted by the Fermi function
        else:
            # Fermi function:
            f = 1/(np.exp(w/T)+1)
            hartree_entries = np.einsum('in,in,n->i',v[:self.n_dof],np.conj(v[:self.n_dof]),f,optimize=True) 
            fock_entries = -np.einsum('in,in,n->i',tmp0,tmp1,f,optimize=True)
            gorkov_entries = (1/2)*np.einsum('in,in,n->i',tmp0,tmp1a,2*f-1,optimize=True)
        return hartree_entries, fock_entries, gorkov_entries

    def anomalous_fourier_transform(self,v):
        coeff=self.coeff
        coeff=np.block([[coeff,np.conj(coeff)],[np.conj(coeff),coeff]])
        return np.dot(coeff.T,v)

    def _set_fields(self, hartree_entries, fock_entries, gorkov_entries):
        """Returns Hartree, Fock, Gorkov"""
        hartree = -np.einsum('ij,j->i',self._hubbard_u,hartree_entries,optimize=True)
        fock =   +np.multiply(self.U_entries,fock_entries)
        gorkov = -np.multiply(self.U_entries,gorkov_entries)

        #free_energy=0
        #free_energy = np.einsum('i,i->',self.U_entries,gorkov_)
        #free_energy+= np.einsum('ij,i,j->',self._hubbard_u,hartree,hartree)
        # free_energy=free_energy+np.sum(self.w[:self.n_dof])
        #print(free_energy)

        return hartree, fock, gorkov
        # return np.real_if_close(hartree), np.real_if_close(fock), np.real_if_close(gorkov)

    def __iter__(self):

        self.iterations=0

        self._hartree_record=[]
        self._fock_record=[]
        self._gorkov_record=[]
        
        return self

    def __next__(self):
        
        hartree, fock, gorkov = self.hartree, self.fock, self.gorkov
        #hartree, fock, gorkov = self._hartree+self.external_hartree, self._fock+self.external_fock, self._gorkov+self.external_gorkov
        
        ########### record fields: ###########
        n_h = len(self._hartree_indices)
        n_f = len(self._fock_indices)
        n_g = len(self._gorkov_indices)
        n_u = np.shape(self.hubbard_indices)[0]
        hubbard_indices=np.transpose(self.hubbard_indices)
        # for i in range(n_u):

        temp=[]
        for index in self._hartree_indices:
            temp.append(complex(self.hartree[index]))
        self._hartree_record.append(temp)
        
        temp_fock=[]
        temp_gorkov=[]

        for i,hubbard_index in enumerate(hubbard_indices):

            for index in self._fock_indices:
                if np.array_equal(index,hubbard_index):
                    temp_fock.append(fock[i])

            for index in self._gorkov_indices:
                if np.array_equal(index,hubbard_index):
                    temp_gorkov.append(gorkov[i])

        self._fock_record.append(temp_fock)
        self._gorkov_record.append(temp_gorkov)

        if np.sum(self._hartree_print_indices)>0:
            print('\nHartree field:')
            temp=np.array(self._hartree_record[-1])
            print(temp[self._hartree_print_indices])
        if np.sum(self._fock_print_indices)>0:
            print('\nFock field:')
            temp=np.array(self._fock_record[-1])
            print(temp[self._fock_print_indices])
        if np.sum(self._gorkov_print_indices)>0:
            print('\nGorkov field:')
            temp=np.array(self._gorkov_record[-1])
            print(temp[self._gorkov_print_indices])


        ##################################

        self.set_mean_field_hamiltonian()
        w,v = self.solve()
        #i=int(self.n_dof/2)
        Eg=np.sum(la.eigh(self._tb_ham,eigvals_only=True))/self.n_dof
        h=np.einsum('ij,i,j->',self._hubbard_u,hartree,hartree,optimize=True)
        f=np.einsum('i,i,i->',self.U_entries,fock,fock,optimize=True)
        g=np.einsum('i,i,i->',self.U_entries,gorkov,gorkov,optimize=True)
        o=np.dot(np.diag(self._tb_ham),hartree)
        n=np.dot(self._tb_ham[self.hubbard_indices],fock)
        #anom=ham[tuple([self.hubbard_indices[0]+self.n_dof,self.hubbard_indices[1]])]
        #gg=np.dot(anom,gorkov)
        if self.iterations>1:
            self.energy_old=np.copy(self.energy)
        self.energy=h-f+g+o-n#+gg
        self.energy=self.energy/self.n_dof
        #print(self.energy)

        if self.iterations>1:
            self.en_old=np.copy(self.en)
        self.en=np.sum(w[self.n_dof:])/self.n_dof
        #print(self.en)

        #if self.iterations>1:
        #    print(self.en<self.en_old)
        if self.temperature==0:
            f_mf=0
        else:
            f_mf =-self.temperature*np.sum(np.log(1+np.exp(-w[self.n_dof:]/self.temperature)))/self.n_dof
        #print(f_mf)

        self.free_energy=self.energy+f_mf
        #print(free_energy)
        #print('------')

        #print(np.sum(self.w[:self.n_dof]))

        hartree_entries, fock_entries, gorkov_entries = self._thermal_density_matrix(w,v)

        hartree, fock, gorkov = self._set_fields(hartree_entries, fock_entries, gorkov_entries)

        hartree = (1-self.friction)*hartree + self.friction*self.hartree
        fock = (1-self.friction)*fock + self.friction*self.fock
        gorkov = (1-self.friction)*gorkov + self.friction*self.gorkov
        
        # update fields:
        self.hartree, self.fock, self.gorkov = hartree, fock, gorkov 
        
        print(np.sum(w[:self.n_dof]))
        return w,v

    def self_consistent_calculation(self):
        """If dos=True, the density of states are calculated once the self consistent
        loop has converged."""

        t = time.time()

        self._set_hubbard_indices()

        self.fock = self.fock[self.hubbard_indices]
        self.gorkov = self.gorkov[self.hubbard_indices]

        iteration = iter(self)
        
        for i in tqdm(range(self.max_iterations)):
            
            self.iterations+=1

            hartree_old = self.hartree
            fock_old = self.fock
            gorkov_old = self.gorkov

            w,v = next(iteration)
            
            eps = self.absolute_convergence_factor
            if (np.allclose(hartree_old,self.hartree,atol=eps) & np.allclose(fock_old,self.fock,atol=eps) & np.allclose(gorkov_old,self.gorkov,atol=eps)):
                break
            if i+1 == self.max_iterations:
                print('Did not converge within max_iterations!')
                break

        self.exec_time = time.time() - t

        self._hartree_record = np.transpose(self._hartree_record)
        self._fock_record = np.transpose(self._fock_record)
        self._gorkov_record = np.transpose(self._gorkov_record)
        
        # Trash non
        del(self._tb_ham)
        del(self._hubbard_u)

        return w,v

    def hartree(self):
        temp=np.reshape(self._hartree,self.extended_dimensions,'F')
        return np.real_if_close(temp)

    def fock(self):
        extended_dimensions=np.append(self.extended_dimensions,self.extended_dimensions)
        temp=np.zeros([self.n_dof,self.n_dof],dtype=complex_dtype)
        temp[self.hubbard_indices]=self._fock
        temp=np.reshape(temp,extended_dimensions,'F')
        return np.real_if_close(temp)

    def gorkov(self):
        extended_dimensions=np.append(self.extended_dimensions,self.extended_dimensions)
        temp=np.zeros([self.n_dof,self.n_dof],dtype=complex_dtype)
        temp[self.hubbard_indices]=self._gorkov
        temp=np.reshape(temp,extended_dimensions,'F')
        return np.real_if_close(temp)

class greens_function(bogoliubov_de_gennes):

    def __init__(self, model: "tight binding or bogoliubov_de_gennes", energy_interval: '1D array', resolution: complex, eigenvalues, eigenvectors, anomalous=False):
        """Calls LAPACK Hermitian matrix solver from scipy.linalg.eigh
        with overwrite to conserve memory

        Attributes
        ----------

        dt : float
            simulation time

        Returns 
        ----------
        dos : array
            density of states
        """
        
        self.extended_dimensions = model.extended_dimensions
        self.n_dof = model.n_dof

        self.energy_interval=energy_interval
        self.resolution=resolution
        
        if anomalous:
            density_matrix = self._anomalous_density_matrix(eigenvectors) 
        else:
            density_matrix = self._density_matrix(eigenvectors) 
        density_of_states = DOS(energy_interval,resolution,eigenvalues,density_matrix)
        n_energy=len(energy_interval)
        self.density_of_states = np.reshape(density_of_states, np.append(self.extended_dimensions,[n_energy]), 'F')

    def local_density_of_states(self, energy, orbital, trace_over_spin=True):
        temp = self.density_of_states
        if trace_over_spin:
            temp = np.sum(temp,-3)
        index = FindNearestValueOfArray(self.energy_interval,energy)
        temp = temp[...,index]
        temp = temp[..., orbital]
        return temp
        # return local_density_of_states(self.density_of_states, index, trace_over)

    def spectrum(self, locations, orbital, trace_over_spin=True):
        temp = self.density_of_states
        if trace_over_spin:
            temp = np.sum(temp,-3)
        temp = temp[..., orbital]
        temp = [temp[location] for location in locations]
        return temp

class Lattice(tight_binding):
    def __init__(self, fig, ax, model, ldos, layer, spins, orbs, annotate_orbs=True, show_cell_borders=True, show_basis_vectors=True, show_hopping=True, show_impurities=True, show_only_centre=True):

        super().__init__(model.dimensions, model.n_spins, model.basis, model.orbitals, model.hoppings, model.impurities)

        self.fig = fig
        self.ax = ax
        # Colors and labels:
        self.color_orb = plt.cm.rainbow(np.linspace(0,1,self.n_orbs))
        self.label_basis = list(string.ascii_uppercase)
        self.label_orb = list(string .ascii_lowercase)

        # Show hopping:
        if show_hopping:
            self.label_hopping=[r'$\delta$']
            for link in self.hoppings:
                self._plt_hopping(*link, show_only_centre)

        # Find orbitals in layer to be plotted:
        for o, orb in enumerate(self.orbitals):
            
            tmp=layer+(spins,o,slice(2))
            temp = self.pos_orb[tmp]
            pos = temp.reshape([self.dimensions[0]*self.dimensions[1]*self.n_spins, 2])
            size = ldos[:,:,:,o].reshape([self.dimensions[0]*self.dimensions[1]*self.n_spins])
            xx, yy = zip(*pos)
            s=size*20/(np.mean(size))
            self.ax.scatter(xx, yy, s=s, color=self.color_orb[o])

        # Impurities
        if show_impurities:
            [V, impurity_loc, impurity_spin, impurity_orb] = self.impurities
            ll=[]
            for loc in impurity_loc:
                if len(loc)==2:
                    loc.append(0)
                ll.append(loc)
            ll = np.array(ll)
            ss=[]
            for spin in impurity_spin:
                if spin in spins:
                    ss.append(spin)
            oo=[]
            for orb in impurity_orb:
                oo.append(orb)
            n_imp = np.shape(ll)[0]
            for o in oo:
                if ll!=[]:
                    pos = (self.pos_orb[ll[:,0],ll[:,1],ll[:,2],ss,o,:2])
                    xx, yy = zip(*pos)
                    s=10/np.mean(self.dimensions[0:1])
                    self.ax.scatter(xx, yy, s=6*s, color='k', marker='o')
                    self.ax.scatter(xx, yy, s=s, color=self.color_orb[o], marker='*')
        
        # Cell borders:
        if show_cell_borders==True:
            width = self.basis[0][0]
            # width = 2*basis[0][0] # 2 for hexagonal self
            height = self.basis[1][1]

            tmp=layer+(0,0,slice(2))
            temp = self.pos_orb[tmp]
            for a_x, a_y in zip(*(temp.reshape([self.dimensions[0]*self.dimensions[1],2])+self.com[:2]).T):
                self.ax.add_patch(Rectangle(
                    xy=(a_x-width/2, a_y-height/2) ,width=width, height=height,
                    linewidth=1, color='blue', fill=False))
        
        # Annotate orbitals:
        for o, orb in enumerate(self.orbitals):
            if annotate_orbs==True:
                self.ax.annotate(self.label_orb[o], self.orbitals[o][:2]+[0.05,0.05], horizontalalignment='right', verticalalignment='bottom', color=self.color_orb[o])
        # Draw basis vectors:
        if show_basis_vectors==True:
            for i in range(2): # Two basis vectors are plotted in the plane
                self.ax.annotate(text='', xytext=self.basis[i][:2],xy=self.orbitals[0][:2],  arrowprops=dict(arrowstyle='<-', lw=2))
                self.ax.annotate(text=self.label_basis[i], xytext=0.5*self.basis[i][:2]+[0.05,0.05],xy=self.orbitals[0][:2])

    def _plt_hopping(self, t: complex, orb_i: int, orb_f: int, cell_hop: tuple, PBC=True, show_only_centre=True):
        if len(cell_hop)==2:
            cell_hop.append(0)
        cell_hop = np.array(cell_hop)
        cells = self.n_cells
        s=0
        if (cell_hop==np.zeros([3])).all() and PBC==True:
            if show_only_centre:
                index=[0]
            else:
                index=np.arange(cells)
            for i in index:
                [xi,yi,zi] = Z_Zn(i, self.dimensions)
                [xf,yf,zf] = np.mod([xi,yi,zi] + cell_hop, self.dimensions)
                orb_0 = self.pos_orb[xi,yi,zi,s,orb_i] 
                orb_1 = self.pos_orb[xf,yf,zf,s,orb_f]
                orb_0=orb_0[:2]
                orb_1=orb_1[:2]
                self.ax.annotate(text='', xytext=orb_1, xy=orb_0,  arrowprops={'arrowstyle': '<->', 'ls': 'dashed'})
                # self.ax.annotate(text=self.label_hopping[0], xytext=0.5*orb_1, xy=orb_0)
        # else:
        for dimension in range(self.n_dim):
            temp=np.zeros(3)
            temp[dimension]=1
            if (temp==cell_hop).all():
                if show_only_centre:
                    index=[0]
                else:
                    index=np.arange(cells)
                for index_i in index: 
                    [xi,yi,zi] = Z_Zn(index_i, self.dimensions)
                    [xf,yf,zf] = np.mod([xi,yi,zi] + cell_hop, self.dimensions)
                    coord_i = self.coord_cell[index_i]
                    edge = bool(coord_i[dimension]>=self.centre[dimension])
                    if edge and not PBC:
                        pass
                    else:
                        orb_0 = self.pos_orb[xi,yi,zi,s,orb_i] 
                        orb_1 = self.pos_orb[xf,yf,zf,s,orb_f]
                        orb_0=orb_0[:2]
                        orb_1=orb_1[:2]
                        self.ax.annotate(text='', xytext=orb_1, xy=orb_0,  arrowprops={'arrowstyle': '<->', 'ls': 'dashed'})

class LocalDensityOfStates(tight_binding):
    def __init__(self, fig, ax, omegas, ldos, omega):

        self.fig = fig
        self.ax = ax

        index = FindNearestValueOfArray(omegas, omega)
        omega = omegas[index]

        eV=np.real(omega)
        epsilon=np.imag(omega)

        self.ldos=ldos

        x,y=np.shape(ldos)[:2]
        
        self.text= (f'DOS map'
        '\n'
        f'$\omega={eV:.2f}$'
        '\n'
        f'$\epsilon={epsilon}$'
        )    
        # k_F = Fermi_vector(mu, t)
        # lambda_F = Friedel_wavelength(k_F)
        # text_fermi = (f'$\lambda_F={lambda_F:.2f}$')
        # text = ('tight_binding parameters: '
               # f'$\mu={mu:.2f}$, '
               # f'$s={s:.3f}$, '
               # f'$\delta={delta:.3f}$, '
               # f'$t={t:.2f}$, '
               # f'$d=({d[0]:.2f},{d[1]:.2f},{d[2]:.2f})$, '
               # f'$N_x = {n_sites_x}$, '
               # f'$N_y = {n_sites_y}$'
               # '\n'
               # 'Impurity location: '
               # f'{impurity_loc}, '
               # f'$V={V:.2f}$'
               # )
        self.interpolation = 'none'
        self.vmin, self.vmax = np.amin(ldos), np.amax(ldos)
        self._extent = [-x/2,x/2,-y/2,y/2]
        self.cmap = ListedColormap(cm.afmhot(np.linspace(0, 0.75, 256)))

    def imshow(self,layer=(ALL,ALL,0)):
        x=self.ldos[layer]
        x=np.fft.fftshift(x)
        self.im = self.ax.imshow(
                x.T, extent=self._extent,
                interpolation=self.interpolation,
                vmin=self.vmin, vmax=self.vmax,
                cmap=self.cmap)

        self.ax.set(xlabel='$x/a_x$', ylabel='$y/a_y$')
        
        self.fig.subplots_adjust(hspace=0.1, wspace=0.1)
        return self.fig,self.ax
    #    cbar.ax.set_title(r'eV')
        
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

    def plot_3D(self):
        import plotly.graph_objects as go
        dimensions=self.dimensions
        X, Y, Z = np.mgrid[0:dimensions[0], 0:dimensions[1], 0:dimensions[2]]
        values = np.fft.fftshift(self.ldos)
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

        return fig

class Energy(tight_binding):
    def __init__(self, fig, ax, model, dos):

        super().__init__(model.dimensions, model.n_spins, model.basis, model.orbitals, model.pbc)
        
        self.fig = fig
        self.ax = ax

        self.dos=dos
        
        self.text= (f'Energy')    

        self.vmin, self.vmax = np.amin(dos), np.amax(dos)

    def plot(self):
        x=np.arange(len(self.dos))
        y=self.dos
        y=np.fft.fftshift(y)
        self.scatter = self.ax.scatter(
                x,y) 
                #vmin=self.vmin, vmax=self.vmax,
                #cmap=self.cmap)

        self.ax.set(xlabel='$x/a_x$', ylabel='Energy')
        
        self.fig.subplots_adjust(hspace=0.1, wspace=0.1)
        return self.fig, self.ax

