#!/usr/bin/env python

# Copyright under GNU General Public License 2021
# by Henry Sheehy

import sys
import os
ROOT_DIR = os.path.dirname(
             os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT_DIR)
sys.path.append('src/')
from imports import *
'''A library for multiorbital non-unitary spin-triplet superconductivity on a lattice
simulation.'''
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

        dim=self.dimensions
        self.edge = np.floor(0.5*np.array(dim),dtype=float)
        for i in range(len(dim)):
            if self.pbc[i]:
                self.edge[i]=+100000000

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

        self.hamiltonian = np.zeros([self.n_dof, self.n_dof], dtype=complex_dtype) 

        self.coeff = self._fourier_coeffs()

    def set_temperature(self,temperature):
        self.temperature=temperature

    #################################### 
    ########### Hamiltonian ############
    ####################################

    def set_onsite_chemical_potential(self, chemical_potential, orbital=None, axes=None):
        # axes is a slice
        spin_matrix = chemical_potential*np.eye(self.n_spins)
        self.set_onsite_spin(spin_matrix, orbital, axes)

    def set_onsite_spin(self, spin_matrix, orbital=None, axes=None):
        # axes is a slice
        if type(orbital)==type(None):
            orbital=ALL
        orb = np.zeros([self.n_orbs])
        orb[orbital] = 1
        orb = np.diag(orb)
        SO_onsite=np.kron(spin_matrix,orb)
        self.set_onsite_spin_orbital(SO_onsite,axes)

    def set_onsite_spin_orbital(self, spin_orbital_onsite, axes=None):
        # axes is a slice
        if type(axes)==type(None):
            axes=[ALL for i in range(self.n_dim)]
        axes=tuple(axes)
        index = Index(self.dimensions)
        index = (index[axes])
        sites = np.zeros([self.n_cells])
        sites[index] = 1
        sites=np.diag(sites)
        self.hamiltonian += np.kron(spin_orbital_onsite,sites)
    
    def set_hopping(self, t: complex, orb_i: int, orb_f: int, cell_hop: tuple, skip_single_site=True):
        spin_tensor=t*np.eye(self.n_spins)
        self.set_spin_hopping(spin_tensor,orb_i,orb_f,cell_hop,skip_single_site)
        
    def set_spin_hopping(self, spin_tensor: '2D-array', orb_i: int, orb_f: int, cell_hop: tuple, skip_single_site=True):
        orbit_tensor=np.zeros([self.n_orbs, self.n_orbs])
        orbit_tensor[orb_i,orb_f]=1
        SO_tensor=kron(spin_tensor,orbit_tensor)
        self.set_spin_orbit_hopping(SO_tensor,cell_hop,skip_single_site)

    def set_spin_orbit_hopping(self, spin_orbit_tensor: '2D-array', cell_hop: tuple, skip_single_site=True):
        self.hamiltonian += self._spin_orbit_hopping(spin_orbit_tensor, cell_hop, skip_single_site)

    def _spin_orbit_hopping(self, spin_orbit_tensor, hop: tuple, skip_single_site=True):
        dim=self.dimensions
        edge=np.array(self.edge,dtype=int)
        pbc=self.pbc
        
        dimension_along_hop = np.dot(dim,np.abs(hop)) #abs for negative direction hopping
        single_site = (dimension_along_hop==1)
        if single_site and skip_single_site:
            temp=np.zeros([self.n_dof,self.n_dof])
            return temp
        temp=np.zeros([self.n_cells,self.n_cells])
        x=np.indices(dim)
        x=np.moveaxis(x,0,-1)
        x=np.add(np.mod(np.add(x, self.centre), dim), -self.centre)
        y=np.copy(x)
        y=y+hop
        # cut out hop outside of boundary:
        indices=np.invert(np.any(y>edge,axis=-1))
        x = x[indices]
        y = y[indices]
        x=coordinates_to_indices(x,dim)
        y=coordinates_to_indices(y,dim)
        x=x.flatten()
        y=y.flatten()
        temp[x,y]+=1
        temp=kron(spin_orbit_tensor,temp)
        temp+=dagger(temp)
        return temp

    def set_impurities(self, V, locations, spin_matrix=None, orbital=None):
        for location in locations:
            self.set_onsite_chemical_potential(V, orbital, location)

    def set_magnetic_impurities(self, M, locations, spin_matrix=None, orbital=None):
        spin_matrix=np.einsum('i,ijk->jk',M,Pauli_vec)
        self.set_onsite_spin(spin_matrix, orbital, locations)

    def unravel_hamiltonian(self):
        self.hamiltonian = np.reshape(self.hamiltonian, np.append(self.extended_dimensions,self.extended_dimensions), 'F')

    def ravel_hamiltonian(self):
        self.hamiltonian = np.reshape(self.hamiltonian, [self.n_dof,self.n_dof], 'F')

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
        if np.shape(self.hamiltonian)==(self.n_dof,self.n_dof):
            ravelled=True
            self.unravel_hamiltonian()
        
        self.hamiltonian = np.fft.ifftn(self.hamiltonian,axes=transform)
        transform = tuple(np.array(transform)+len(self.extended_dimensions))
        self.hamiltonian = np.fft.fftn(self.hamiltonian,axes=transform)
    
        self.hamiltonian = np.fft.fftn(self.hamiltonian,axes=inverse_transform)
        inverse_transform = tuple(np.array(inverse_transform)+len(self.extended_dimensions))
        self.hamiltonian = np.fft.ifftn(self.hamiltonian,axes=inverse_transform)

        if ravelled:
            self.ravel_hamiltonian()
    
        ########################################
        ######### Statistical Mechanics ########
        ########################################

    def _density_matrix(self, eigenvectors):
        return np.einsum('in,in->in',eigenvectors[:self.n_dof],np.conj(eigenvectors[:self.n_dof]),optimize=True)
    
    def solve_eigensystem(self):
        w,v = la.eigh(self.hamiltonian, overwrite_a=True)
        #from scipy.sparse.linalg import eigsh
        # w,v = eigsh(ham, k=50, which='SM', return_eigenvectors=True)
        # v=self.fourier_transform_psi(v)
        # v=self.inv_fourier_transform_psi(v)
        return w,v

    def density_of_states(self, energy_interval: '1D array', resolution: complex, eigenvalues=None, eigenvectors=None):
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
        self.energy_interval=energy_interval
        t = time.time()
        if type(eigenvalues)==type(None) and type(eigenvectors)==type(None):
            # eigenvalues and eigenvectors don't exist... solving
            eigenvalues,eigenvectors=self.solve_eigensystem()
        density_matrix = self._density_matrix(eigenvectors) 
        density_of_states = DOS(energy_interval,resolution,eigenvalues,density_matrix)
        self.exec_time = time.time() - t
        n_energy=len(energy_interval)
        self.density_of_states = np.reshape(density_of_states, np.append(self.extended_dimensions,[n_energy]), 'F')
        return self.density_of_states

    def local_density_of_states(self, density_of_states, energy, trace_over=True):
        index = FindNearestValueOfArray(self.energy_interval,energy)
        return local_density_of_states(density_of_states, index, trace_over)
    
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
        
        self._hubbard_u = None
        self._hubbard_indices = None

        #self._hartree = np.zeros([self.n_dof], dtype=complex_dtype)
        #self._fock = np.diag(np.zeros([self.n_dof], dtype=complex_dtype))
        #self._gorkov = np.diag(np.zeros([self.n_dof], dtype=complex_dtype))

        #self.external_hartree = 0*np.copy(initial_hartree)
        #self.external_fock = 0*np.copy(initial_fock)
        #self.external_gorkov = 0*np.copy(initial_gorkov)

        #self.external_hartree = external_hartree
        #self.external_fock = external_fock
        #self.external_gorkov = external_gorkov


        self._hartree_indices = []
        self._fock_indices = []
        self._gorkov_indices = []
        
    def set_hamiltonian(self,hartree,fock,gorkov):
        
        n_dof = self.n_dof

        hamiltonian = np.zeros([2*n_dof,2*n_dof],dtype=complex_dtype)
        
        hamiltonian[:n_dof,:n_dof]=self._tb_ham+np.diag(hartree)
        hamiltonian[n_dof:,n_dof:]=-np.conj(self._tb_ham)-np.conj(np.diag(hartree))
        hamiltonian[self.hubbard_indices[0],self.anomalous_indices[1]]=fock
        hamiltonian[self.anomalous_indices[0],self.hubbard_indices[1]]=-np.conj(fock)
        hamiltonian[self.hubbard_indices[0],self.anomalous_indices[1]]=-np.conj(gorkov)
        hamiltonian[self.anomalous_indices[0],self.hubbard_indices[1]]=gorkov

        return hamiltonian

    def set_fields(self,hartree=None,fock=None,gorkov=None):
        self._set_hubbard_indices()
        self.hartree = hartree.astype(complex_dtype)
        self.fock = fock[self.hubbard_indices].astype(complex_dtype)
        self.gorkov = gorkov[self.hubbard_indices].astype(complex_dtype)

    def set_external_fields(self,external_hartree=None,external_fock=None,external_gorkov=None):
        self.external_hartree = external_hartree
        self.external_fock = external_fock
        self.external_gorkov = external_gorkov

    def set_hubbard_U(self, spin_matrix: 'array', orb_i: int, orb_f: int, cell_hop: tuple, skip_single_site=True):
        orbit_tensor=np.zeros([self.n_orbs, self.n_orbs])
        orbit_tensor[orb_i,orb_f]=1
        SO_tensor=kron(spin_tensor,orbit_tensor)
        self.set_hubbard_u_spin_orbit(SO_tensor,cell_hop,skip_single_site)
        
    def set_hubbard_U_spin_orbit(self, spin_orbit_tensor: '2D-array', cell_hop: tuple, skip_single_site=True):
        if type(self._hubbard_u)==type(None):
            self._hubbard_u=np.zeros([self.n_dof, self.n_dof], dtype=complex_dtype)
        self._hubbard_u += self._spin_orbit_hopping(spin_orbit_tensor, cell_hop, skip_single_site)

    def _set_hubbard_indices(self):

        self.hubbard_indices = np.nonzero(self._hubbard_u)
        self.U_entries = self._hubbard_u[self.hubbard_indices]

        self.anomalous_indices=np.array(self.hubbard_indices)
        self.anomalous_indices+=self.n_dof
        self.anomalous_indices=tuple(self.anomalous_indices)

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

    def anomalous_density_of_states(self, energy_interval: '1D array', resolution: complex, eigenvalues=None, eigenvectors=None):

        self.energy_interval=energy_interval
        t = time.time()
        if type(eigenvalues)==type(None) and type(eigenvectors)==type(None):
            # eigenvalues and eigenvectors don't exist... solving
            eigenvalues,eigenvectors=self.solve_eigensystem()
        anomalous_density_matrix = self._anomalous_density_matrix(eigenvectors) 
        density_of_states = DOS(energy_interval,resolution,eigenvalues,anomalous_density_matrix)
        self.exec_time = time.time() - t
        n_energy=len(energy_interval)
        self.density_of_states = np.reshape(density_of_states, np.append(self.extended_dimensions,[n_energy]), 'F')
        return self.anomalous_density_of_states

    def anomalous_fourier_transform(self,v):
        coeff=self.coeff
        coeff=np.block([[coeff,np.conj(coeff)],[np.conj(coeff),coeff]])
        return np.dot(coeff.T,v)

    def solve(self, hamiltonian):
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

        t = time.time()

        w,v = la.eigh(hamiltonian, overwrite_a=True)

        self.dos, self.ados = self.density_of_states(w,v), self.anomalous_density_of_states(w,v)
        self.exec_time = time.time() - t

        return self.dos, self.ados


    def _set_fields(self, hartree_entries, fock_entries, gorkov_entries):
        """Returns Hartree, Fock, Gorkov"""
        hartree = -np.einsum('ij,j->i',self._hubbard_u,hartree_entries,optimize=True)
        fock =   +np.multiply(self.U_entries,fock_entries)
        gorkov = -np.multiply(self.U_entries,gorkov_entries)
        print(hartree[0])
        print(fock[0])
        print(gorkov[0])

        #free_energy=0
        #free_energy = np.einsum('i,i->',self.U_entries,gorkov_)
        #free_energy+= np.einsum('ij,i,j->',self._hubbard_u,hartree,hartree)
        # free_energy=free_energy+np.sum(self.w[:self.n_dof])
        #print(free_energy)

        return hartree, fock, gorkov
        # return np.real_if_close(hartree), np.real_if_close(fock), np.real_if_close(gorkov)

    def set_field_tracking(self, hartree_indices=[], fock_indices=[], gorkov_indices=[]):
        temp=[]
        for index in hartree_indices:
            temp.append(index)    
        self._hartree_indices=temp
        temp=[]
        for index in fock_indices:
            tmp=[]
            tmp.append(np.argwhere(np.logical_and((self.hubbard_indices)[0]==index[0], (self.hubbard_indices)[1]==index[1])))
            tmp=np.concatenate(tmp).ravel().tolist()[0]
            temp.append(tmp)
        self._fock_indices=temp
        temp=[]
        for index in gorkov_indices:    
            tmp=[]
            tmp.append(np.argwhere(np.logical_and((self.hubbard_indices)[0]==index[0], (self.hubbard_indices)[1]==index[1])))
            tmp=np.concatenate(tmp).ravel().tolist()[0]
            temp.append(tmp)
        self._gorkov_indices=temp

    def __iter__(self):

        self.iterations=0

        self._hartree_list=[]
        self._fock_list=[]
        self._gorkov_list=[]
        
        return self

    def __next__(self):
        
        hartree, fock, gorkov = self.hartree, self.fock, self.gorkov
        #hartree, fock, gorkov = self._hartree+self.external_hartree, self._fock+self.external_fock, self._gorkov+self.external_gorkov
        
        ########### tracking ###########
        temp=[]
        for index in self._hartree_indices:
            temp.append(hartree[index])
        self._hartree_list.append(temp)

        temp=[]
        for index in self._fock_indices:
            temp.append(fock[index].tolist())
        self._fock_list.append(temp)

        temp=[]
        for index in self._gorkov_indices:
            temp.append(gorkov[index].tolist())
        self._gorkov_list.append(temp)
        ##################################
        ham = self.set_hamiltonian(hartree,fock,gorkov)
        w,v = la.eigh(ham, overwrite_a=True)
        #i=int(self.n_dof/2)

        Eg=np.sum(la.eigh(self._tb_ham,eigvals_only=True))/self.n_dof
        h=np.einsum('ij,i,j->',self._hubbard_u,hartree,hartree,optimize=True)
        f=np.einsum('i,i,i->',self.U_entries,fock,fock,optimize=True)
        g=np.einsum('i,i,i->',self.U_entries,gorkov,gorkov,optimize=True)
        o=np.dot(np.diag(self._tb_ham),hartree)
        n=np.dot(self._tb_ham[self.hubbard_indices],fock)
        anom=ham[tuple([self.hubbard_indices[0]+self.n_dof,self.hubbard_indices[1]])]
        gg=np.dot(anom,gorkov)
        if self.iterations>1:
            self.energy_old=np.copy(self.energy)
        self.energy=h-f+g+o-n+gg
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

    def self_consistent_calculation(self, energy_interval=None, resolution=None):
        """If dos=True, the density of states are calculated once the self consistent
        loop has converged."""

        t = time.time()

        iteration = iter(self)
        
        self._tb_ham = self.hamiltonian

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

        if type(energy_interval)!=type(None) and type(resolution)!=type(None):
            # v=self.anomalous_fourier_transform(v)
            self.dos, self.ados = self.density_of_states(energy_interval,resolution,w,v), self.anomalous_density_of_states(energy_interval,resolution,w,v)
        self.exec_time = time.time() - t

        self._hartree_list = np.transpose(self._hartree_list)
        self._fock_list = np.transpose(self._fock_list)
        self._gorkov_list = np.transpose(self._gorkov_list)
        
        # Trash non
        del(self._tb_ham)
        del(self._hubbard_u)

        return self.hartree, self.fock, self.gorkov, self.iterations

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

     
    # if Include_DOS==False:
    #     with open(out_folder+confname+'.npz', 'wb') as f:
    #         np.savez(f,Hartree=Hartree,Fock=Fock,Gorkov=Gorkov,iterations=iterations,executionTime=executionTime)
    #         print(iterations)
    # else:
    # print('Running DOS')
    # density_matrix = lattice._density_matrix(v)
    # dos = DOS(omegas, w, density_matrix)
    # density_matrix = lattice._anomalous_density_matrix(v)
    # ados = DOS(omegas, w, density_matrix)
    # with open(confname+'.npz', 'wb') as f:
    #     np.savez(f,
    #     Hartree=Hartree,Fock=Fock,Gorkov=Gorkov,
    #     iterations=iterations,
    #     executionTime=executionTime,
    #     dos=dos,
    #     ados=ados)

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
    def __init__(self, fig, ax, model, ldos, omega):

        super().__init__(model.dimensions, model.n_spins, model.basis, model.orbitals, model.hoppings, model.impurities)
        
        self.fig = fig
        self.ax = ax

        omegas = model.omegas

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
                x, extent=self._extent,
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

