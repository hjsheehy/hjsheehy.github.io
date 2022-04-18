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
    return np.reshape(np.arange(np.prod(dimensions)),dimensions,'C')

def coordinates_to_indices(array, dimensions):
    dim=len(dimensions)
    temp=np.mod(array,dimensions)
    multipler=np.ones([dim],dtype=int)
    for i in range(1,dim):
        multipler[i:]*=dimensions[i-1]
    temp=np.multiply(temp,multipler)
    temp=np.sum(temp,axis=-1)
    return temp

def xpts(dimensions):
    dim=len(dimensions)
    tmp=np.indices(dimensions,dtype=float)
    for i in range(dim):
        tmp[i]=tmp[i]-int(dimensions[i]/2)
        tmp[i]=np.roll(tmp[i],-int(dimensions[i]/2),axis=i)
    temp=[]
    for i in range(dim):
        x=tmp[i].flatten()
        temp.append(x)
    tmp=np.array(temp)
    return tmp

def kpts(dimensions):
    dim=len(dimensions)
    tmp=np.indices(dimensions,dtype=float)
    for i in range(dim):
        tmp[i]=tmp[i]-int(dimensions[i]/2)
        tmp[i]=np.roll(tmp[i],-int(dimensions[i]/2),axis=i)
    temp=[]
    for i in range(dim):
        x=tmp[i].flatten()
        x=x/(dimensions[i])
        temp.append(x)
    tmp=np.array(temp)
    tmp=2*np.pi*tmp
    return tmp

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

###################################################################
########################### Orbitals ##############################
###################################################################
class Orbital():
    """An atomic orbital, instantiated with a position and name"""

    def __init__(self,name:str,label=None):
        self.name=name
        self.label=label
        self.Z=1
        self.index=None

        if name=='1s':
            n=1
            label='$1s$'

            def radial_wave_function(self,r:float)->float:
                return 2*self.Z**(3/2)*np.exp(-self._rho(r)/2)

            def angular_wave_function(self)->float:
                return np.sqrt(1/(4*np.pi))

            def wavefunction(self,x:float,y:float,z:float)->float:
                r=self.radius(x,y,z)
                return self.radial_wave_function(r)*self.angular_wave_function()

        if name=='d_xz':
            n=4
            label='$d_{xz}$'

            def radial_wave_function(self,r:float)->float:
                rho=self._rho(r)
                return (1/(96*np.sqrt(5)))*(6-rho)*rho**2*Z**(3/2)*np.exp(-rho/2)

            def angular_wave_function(self,r:float,x:float,z:float)->float:
                return np.sqrt(60/(16*np.pi))*x*z/r**2

            def wavefunction(self,x:float,y:float,z:float)->float:
                r=self.radius(x,y,z)
                return self.radial_wave_function(r)*self.angular_wave_function(r,x,z)

        if name=='d_yz':
            n=4
            label='$4d_{yz}$'

            def radial_wave_function(self,r:float)->float:
                rho=self._rho(r)
                return (1/(96*np.sqrt(5)))*(6-rho)*rho**2*z**(3/2)*np.exp(-rho/2)

            def angular_wave_function(self,r:float,y:float,z:float)->float:
                return np.sqrt(60/(16*np.pi))*y*z/r**2

            def wavefunction(self,x:float,y:float,z:float)->float:
                r=self.radius(x,y,z)
                return self.radial_wave_function(r)*self.angular_wave_function(r,y,z)
          
        if name=='d_xy':
            n=4
            label='$4d_{xy}$'

            def radial_wave_function(self,r:float)->float:
                rho=self._rho(r)
                return (1/(96*np.sqrt(5)))*(6-rho)*rho**2*Z**(3/2)*np.exp(-rho/2)

            def angular_wave_function(self,r:float,x:float,y:float)->float:
                return np.sqrt(60/(16*np.pi))*x*y/r**2

            def wavefunction(self,x:float,y:float,z:float) -> float:
                r=self.radius(x,y,z)
                return self.radial_wave_function(r)*self.angular_wave_function(r,x,y)

    def __repr__(self):
        return "Orbital()"

    def __str__(self):
        return "{self.name} orbital"
    
    def set_effective_nuclear_charge(self,Z:float)->float:
        self.Z=Z

    def _rho(self,r:float)->float:
        return 2*self.Z*r/self.n

    def radius(self,x:float,y:float,z:float)->float:
        return numpy.linalg.norm([x,y,z])

    def electron_density(self,x:float,y:float,z:float)->float:
        return np.abs(self.wavefunction(x,y,z))**2

###################################################################
####################### Crystal lattice ###########################
###################################################################
class Atom():
    """An atom, instantiated with a Cartesian position, a name and an empty list of orbitals"""
#     Vectors=list[list[float]]
    def __init__(self, position,name:str):
        self.position=position
        self.name=name
        self._index=None
        self._orbitals=[]
        self._orbital_dict=dict()
        self._orbital_indices=[]
        self._counter=0

    def __repr__(self):
        return "Atom()"

    def __str__(self):
        return f"{self.name} atom, located at {self.position} with {self.n_orbitals} orbitals: {self.orbitals}"
    
    def add_orbital(self, orbital):
        if type(orbital)==str:
            orbital=Orbital(orbital)
        orbital._index=self._counter
        # if not isinstance(orbital,Orbital):
        #     raise ValueError(f'{orbital} is not an instance of the Orbital class!')
        self._orbitals.append(orbital)
        self._orbital_dict[orbital.name]=self._counter
        self._orbital_indices.append(self._counter)
        self._counter+=1

    def _orbital_index(self,orbital):
        if type(orbital)==str:
            return self._orbital_dict[orbital]
        elif type(orbital)==float:
            orbital=int(orbital)
        elif type(orbital)==int or np.issubdtype(orbital,np.integer):
            return orbital%self.n_orbitals
        elif isinstance(orbital,Orbital):
            return self.orbital._index
        else:
            raise ValueError(f'{orbital} not str, int, or Orbital class!')

    def _orbital_name(self,orbital):
        if type(orbital)==float:
            orbital=int(orbital)
        if type(orbital)==int:
            return self.orbitals[orbital%self.n_orbitals]
        elif type(orbital)==str:
            return orbital
        elif isinstance(orbital,Orbital):
            return self.orbital.name
        else:
            raise ValueError(f'{orbital} not str, int, or Orbital class!')

    @property
    def orbitals(self):
        return list(self._orbital_dict.keys())

    @property
    def n_orbitals(self):
        return len(self._orbitals)

    @property
    def _indices(self):
        if type(self._index)!=type(None):
            return np.array(self._orbital_indices)+self._index
        else:return self._orbital_indices

class CrystalLattice():
    """Crystal lattice"""

    def __init__(self,lattice_vectors=[[1,0],[0,1]],name='lattice'):
        self.lattice_vectors=np.array(lattice_vectors)
        self.name=name
        self._atoms=[]
        self._atom_dict=dict()
        self._atom_indices=[]
        self._counter=0
        self.n_spins=2
        n_dimensions=len(self.lattice_vectors)
        self.bulk_calculation=False
        # if np.sum(np.abs(self.lattice_vectors))==0:
        #     self.lattice_vectors=np.array([])
        # if list(self.lattice_vectors)==list([]):
        #     self.lattice_vectors=[]
        # # elif n_dimensions<len(self.lattice_vectors[0]):
        # else:
        #     self.nonzero_axes=np.array(np.nonzero(self.lattice_vectors)).T
        #     self.lattice_vectors=[]
        #     for i in range(len(self.nonzero_axes)):
        #         j=self.nonzero_axes[i][1]
        #         self.lattice_vectors.append([self.lattice_vectors[d,i] for d in range(n_dimensions)])
        # # else:
        # #     self.lattice_vectors=self.lattice_vectors
        #         self.lattice_vectors=np.array(self.lattice_vectors)

        # null_space=la.null_space(self._A)
        # self.lattice_vectors=np.copy(self.lattice_vectors)
        # self.lattice_vectors[null_space]
        # self._A=np.vstack(self.lattice_vectors).T
        # if np.linalg.det(self._A)==0:
        #     raise ValueError('Overcomplete basis!')
        # self._A_inv=np.linalg.inv(self._A)

        self._bulk=[True for nn in range(self.n_dimensions)]
        self._pieces=[1 for nn in range(self.n_dimensions)]
        self._glue_edgs=[True for nn in range(self.n_dimensions)]

    def __repr__(self):
        return "Crystal_lattice()"

    def __str__(self):
        return f"{self.name} lattice, with {self.n_atoms} atoms."

    def cartesian_coordinates(self,fractional_coordinates):
        Y = np.vstack(fractional_coordinates).T
        return np.matmul(self._A, Y.T).T

    def fractional_coordinates(self,cartesian_coordinates):
        Y = np.vstack(cartesian_coordinates).T
        return np.matmul(self._A_inv, Y.T).T

    def add_atom(self, atom):
        if not isinstance(atom,Atom):
            raise ValueError(f'{atom} is not an instance of the Atom class!')
        if atom.name in self.atoms:
            raise SyntaxError(f'Atom name {atom.name} already used. Please use a different one.')
        if len(atom.position)<self.n_dimensions:
            raise ValueError(f'Position of atom {atom} not {self.n_dimensions}-dimensional!')
        atom._index=self._counter
        self._atoms.append(atom)
        self._atom_dict[atom.name]=self._counter
        self._atom_indices.append(self._counter)
        self._counter+=1

    def _atom_index(self,atom):
        if type(atom)==float:
            atom=int(atom)
        if type(atom)==str:
            return self._atom_dict[atom]
        elif type(atom)==int or isinstance(atom,np.int64):
            if atom==0:
                return 0
            return atom%self.n_atoms
        elif isinstance(atom,Atom):
            return atom._index
        else:
            raise ValueError(f'{atom} not str, int, or Atom class!')

    def _atom_name(self,atom):
        if type(atom)==float:
            atom=int(atom)
        if type(atom)==int:
            return self.atoms[atom%self.n_atoms]
        elif type(atom)==str:
            return atom
        elif isinstance(atom,Atom):
            return atom.name
        else:
            raise ValueError(f'{atom} not str, int, or Atom class!')

    def atom(self,atom):
        if type(atom)==float:
            atom=int(atom)
        if isinstance(atom,Atom):
            return atom
        elif isinstance(atom, (np.floating, float, int, np.int_)):
            if atom==0:
                return self._atoms[0]
            return self._atoms[atom%self.n_atoms]
        elif type(atom)==str:
            return self.atom(self._atom_index(atom))
        else:
            raise ValueError(f'{atom} not str, int, or Atom class!')

    def index(self,atom,orbital,spin=0):
        spin=self.spin(spin)
        a=self._atom_index(atom)
        n=self._atoms[a].n_orbitals
        o=self.atom(atom)._orbital_index(orbital)
        o=self._atoms[a]._orbital_index(o)+(spin%self.n_spins)*int(self.n_atoms_orbitals)
        if a==0:
            return o
        else:
            return self.index(a-1,-1)+o+1
    
    @property
    def spins(self):
        if self.n_spins==1:
            return None
        elif self.n_spins==2:
            return { 'up': 0, 'dn': 1, 'down': 1 }
        else:
            raise ValueError(f'{self.n_spins} is nither 1 or 2!')
            
    def spin(self,spin):
        if type(self.n_spins)==int:
            return spin
        elif type(spin)==str:
            return self.spins[spin]
        else:
            raise ValueError(f'{self.n_spins} is neither 1 or 2!')

    def n_orbitals(self,atom):
        return self.atom(atom).n_orbitals
        # return [atom.n_orbitals for atom in self._atoms]

    @property
    def atoms(self):
        return list(self._atom_dict.keys())

    @property
    def atom_positions(self):
        return [atom.position for atom in self._atoms]

    @property
    def n_dimensions(self):
        return len(self.lattice_vectors)
    
    @property
    def n_atoms(self):
        return len(self._atoms)

    @property
    def n_cells(self):
        return np.prod(self._pieces)
    
    @property
    def n_atoms_orbitals(self):
        return np.sum([atom.n_orbitals for atom in self._atoms])

    @property
    def n_atoms_orbitals_spins(self):
        return np.sum([atom.n_orbitals for atom in self._atoms])*self.n_spins

    @property
    def n_dof(self):
        return self.n_spins*self.n_atoms_orbitals*self.n_cells

    @property
    def _extended_dimensions(self):
        return self._pieces+[self.n_atoms_orbitals,self.n_spins]

    @property
    def centre(self):
        """centred-coordinates"""
        return np.array(np.array(self._pieces)/2, dtype=int)  

    @property
    def com(self):
        return np.average(self.atom_positions,axis=0)

    @property
    def edge(self):
        edge = np.floor(0.5*np.array(self._pieces),dtype=float)
        for i in range(self.n_dimensions):
            if self._glue_edgs[i]:
                edge[i]=edge[i]+100000000000
            elif self._pieces[i]%2==0:
                edge-=1
        return edge

    def cut(self, n_cells, axes, glue_edgs=True):
        if type(axes)==int:
            axes=[axes]
        for axis in axes:
            if axis>=self.n_dimensions:
                raise ValueError(f'axis {axis} is greater than dimensions {self.n_dimensions}!')
            self._bulk[axis]=False
            self._pieces[axis]=n_cells
            self._glue_edgs[axis]=glue_edgs

class StatisticalMechanics():
    def __init__(self, eigenvectors, eigenvalues, temperature=0):
        self.eigenvectors=eigenvectors
        self.eigenvalues=eigenvalues
        self.temperature=temperature

    @property
    def density_matrix(self):
        if self.bulk_calculation:
            return np.einsum('kin,kin->kin',self.eigenvectors[:,:self.n_dof],np.conj(self.eigenvectors[:,:self.n_dof]),optimize=True)
        else:
            return np.einsum('in,in->in',self.eigenvectors[:self.n_dof],np.conj(self.eigenvectors[:self.n_dof]),optimize=True)

    @property
    def thermal_density_matrix(self):
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
        if self.n_dof==len(self.eigenvalues):
            return ValueError('_thermal_density_matrix only callable for Bogoliubov de Gennes class, not TightBinding!')

        w,v=self.eigenvalues,self.eigenvectors
        T=self.temperature

        tmp0=v[self._hubbard_indices[0]]
        tmp1=np.conj(v)[self._hubbard_indices[1]]
        tmp1a=v[self._anomalous_indices[1]]
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

class Tightbinding(CrystalLattice, StatisticalMechanics):
    
    def __init__(self,lattice_vectors,name):

        self._model='tb'

        self.kpts = None

        self.temperature=0

        self.hoppings=[]
        self.hopping_labels=[]
        self.impurities=[]
        self.impurity_labels=[]

        self._hamiltonian = None
        self._dos = None
        self._del_eigsys = True

        self._flattened_kpts = None

        super().__init__(lattice_vectors,name)

        self._k_axes = [False for i in range(self.n_dimensions)]

    @property
    def n_kpts(self):
        return list(np.shape(self.kpts[0]))

    @property
    def n_total_kpts(self):
        return np.prod(self.n_kpts)

    @property
    def flattened_kpts(self):
        if type(self._flattened_kpts)==type(None):
            self._flattened_kpts = np.reshape(self.kpts,[self.n_dimensions,self.n_total_kpts])
        return self._flattened_kpts

    def set_temperature(self,temperature):
        self.temperature=temperature

    def _onsite_tensor(self, onsite, atom=None, orbital=None, spin=None, position_coordinates=None):
        """An onsite tensor, e.g. for chemical potential, impurities, spin-flip. 
Specify either the Atom, or the Crystal_lattice. If the Atom is None, then all Atoms in the Crystal_lattice are inputted. If the Orbital is None, then all Orbitals in each Atom are inputted.
The onsite is input as a scalar, a pair (for each spin), a 2-matrix (spin-flips) or 2*orbital-matrix (spin-orbit) i.e. kron(spin,orbit)"""
        indices=[]
        if type(atom)!=type(None):
            _atoms=[self._atom_index(atom)]
        else:
            _atoms=self._atom_indices
        for atom in _atoms:
            n_orbitals=self.n_orbitals(atom)
            if type(orbital)==type(None):
                for orb in np.arange(self.n_orbitals(atom)):
                    indices.append(self.index(atom,orb))
            else:
                indices.append(self.index(atom,orbital))

        if np.isscalar(onsite):
            if self.n_spins==1:
                spin_tensor=onsite
            elif self.n_spins==2:
                if type(spin)==type(None):
                    spin_tensor=onsite*np.eye(2)
                elif type(spin)!=int:
                    raise ValueError(f'{spin} is nither 1 or 2!')
                else:
                    spin_tensor=np.zeros(2)
                    spin_tensor[spin]=onsite
                    spin_tensor=np.diag(spin_tensor)
                    
        elif np.shape(onsite)==(2,):
                spin_tensor=np.diag(onsite)

        elif np.shape(onsite)==(2,2) and self.n_spins==2:
                spin_tensor=onsite

        so=self.n_spins*self.n_orbitals(0)
        if np.array([self.n_orbitals(0)==self.n_orbitals(o) for o in range(self.n_atoms)],bool).all() and np.shape(onsite)==(so,so):
            spin_orbit=onsite
            atom=np.zeros(self.n_atoms)
            atom[self.atoms]=1
            atom=np.diag(atom)
            so_tensor=np.kron(spin_orbit,atom)

        elif np.shape(onsite)==(self.n_atoms_orbitals_spins,self.n_atoms_orbitals_spins):
            so_tensor=onsite

        elif type(atom)!=type(None) and np.shape(onsite)==(self.atom(atom).n_orbitals*self.n_spins,self.atom(atom).n_orbitals*self.n_spins):
            so_tensor=np.zeros([self.n_atoms_orbitals_spins,self.n_atoms_orbitals_spins],dtype=COMPLEX)
            for i in range(len(spins_i)):
                for j in range(len(orbitals_i)):
                    spin1=spins_i[i]
                    orbital1=orbitals_i[j]
                    spin2=spins_f[i]
                    orbital2=orbitals_f[j]
                    index1.append(self.index(atom_i,orbital1,spin1))
                    index2.append(self.index(atom_f,orbital2,spin2))
            so_tensor[index1,index2]=onsite.flatten()

        else:
            so_tensor=np.zeros(self.n_atoms_orbitals,dtype=COMPLEX)
            so_tensor[indices]=1
            so_tensor=np.diag(so_tensor)
            so_tensor=np.kron(spin_tensor,so_tensor)

        if type(position_coordinates)==type(None):
            sites=np.eye(self.n_cells)
            return kron(so_tensor,sites)
        if len(np.shape(position_coordinates))==1:
            position_coordinates=[position_coordinates]
        if len(np.shape(position_coordinates))==2:
            sites=np.zeros([self.n_cells,self.n_cells])
            for x in position_coordinates:
                x=np.moveaxis(x,0,-1)
                x=np.add(np.mod(np.add(x, self.centre), self._pieces), -self.centre)
                x=coordinates_to_indices(x,self._pieces)
                x=x.flatten()
                sites[x,x]=1
            return kron(so_tensor,sites)
        else:
            raise ValueError(f'{position_coordinates} not a list of coordinates!')
        # return kron(sites,so_tensor)

    def _cutout_hopping(self, hop_vector=None):

        if type(hop_vector)==type(None):
            hop_vector=np.zeros(self.n_dimensions)

        hop_vector=np.array(hop_vector)

        if np.logical_and(np.greater_equal(hop_vector,self._pieces),np.invert(self._bulk)).any():
            raise ValueError(f'hop_vector {hop_vector} hops outside the cut!')            

        if np.array_equal(hop_vector,np.zeros([self.n_dimensions])):
            temp = np.eye(self.n_cells)
        else:
            temp=np.zeros([self.n_cells,self.n_cells])
            x=np.indices(self._pieces)
            x=np.moveaxis(x,0,-1)
            x=np.add(np.mod(np.add(x, self.centre), self._pieces), -self.centre)
            y=np.copy(x)
            y=y+hop_vector
            # cut out hop outside of boundary (the edge has been set to finite infinity for pbc):
            indices1=np.invert(np.any(y>self.edge,axis=-1)) # cuts excess to the right
            indices2=np.invert(np.any(y<-self.edge,axis=-1)) # cuts to the left
            indices=np.logical_and(indices1,indices2) # combines
            x = x[indices]
            y = y[indices]
            x=coordinates_to_indices(x,self._pieces)
            y=coordinates_to_indices(y,self._pieces)
            x=x.flatten()
            y=y.flatten()
            temp[x,y]+=1
        return temp

    def _hopping_tensor(self, hopping_amplitude, k=None, atom_i=None, atom_f=None, orbital_i=None, orbital_f=None, hop_vector=None, spin_i=None, spin_f=None, add_time_reversal=True):
        """If TR(hopping)==+\-hopping, then TR not added."""
        if type(k)==type(None):
            temp=self._cutout_hopping(hop_vector)
        else:
            temp=[1]
            hopping_amplitude = hopping_amplitude(k)
        if isinstance(hopping_amplitude, (int, float)):
            # scalar->spin pair
            hopping_amplitude=np.array([hopping_amplitude for i in range(self.n_spins)])
        if np.shape(hopping_amplitude)==(self.n_spins) or np.shape(hopping_amplitude)==(self.n_spins,):
            # spin pair->spin matrix
            hopping_amplitude=np.diag(hopping_amplitude)
        if np.shape(hopping_amplitude)==(self.n_spins,self.n_spins):
            # spin matrix->spin-atom-orbit matrix
            if type(orbital_f)!=type(orbital_f):
                raise ValueError('orbital_i and orbital_f not same type!')
            if type(orbital_i)==type(None):
                if type(atom_i)==type(None):
                    indices_i=np.arange(self.n_atoms_orbitals)
                    indices_f=np.arange(self.n_atoms_orbitals)
                else:
                    indices_i=self.atom(atom_i)._indices
                    indices_f=self.atom(atom_f)._indices
            else:
                indices_i=[self.index(atom_i,orbital_i)]
                indices_f=[self.index(atom_f,orbital_f)]
            
            orbital_matrix=np.zeros(self.n_atoms_orbitals)
            orbital_matrix=np.diag(orbital_matrix)
            orbital_matrix[indices_i,indices_f]=1
            spin_matrix=hopping_amplitude
            hopping_amplitude=np.kron(spin_matrix,orbital_matrix)

        if type(atom_i)!=type(None) and atom_i==atom_f:
            # spin-atom-orbit matrix
            so_dof=self.n_spins*self.atom(atom_i).n_orbitals
            if np.shape(hopping_amplitude)==(so_dof,so_dof):
                # spin-orbit matrix
                A=np.zeros(self.n_atoms_orbitals_spins)
                A=np.diag(A)
                a=self.atom(atom_i)
                ai=a._index
                ao=ai+a.n_orbitals+self.n_spins
                A[ai:ao,ai:ao]=hopping_amplitude
                hopping_amplitude=A

        if np.shape(hopping_amplitude)==(self.n_atoms_orbitals_spins,self.n_atoms_orbitals_spins):
            # unit cell matrix
            pass

        temp=kron(hopping_amplitude,temp)
        # temp=kron(temp,hopping_amplitude)

        # onsite = bool(atom_i==atom_f and (np.array_equal(hop_vector, np.zeros([self.n_dimensions])) or (type(k) is not type(None))))
        TRS=dagger(temp)
        onsite=False
        if (TRS==temp).all() or (TRS==-temp).all():
            onsite=True
        if add_time_reversal and not onsite:
            # temp+=dagger(temp)
            temp+=TRS
        return temp

    def _bulk_hopping(self, k, hop_vector):
        if np.array_equal(hop_vector, np.zeros([self.n_dimensions])):
                return 2
        else:
            return 2*np.cos(np.dot(k,hop_vector))

#     def _bulk_momentum(self, k, func_k, atom_i=None, atom_f=None, orbital_i=None, orbital_f=None, spin_i=None, spin_f=None, add_time_reversal=True):
#         hopping_amplitude = func_k(k)
#         hop_vector = None
#         return self._hopping_tensor(hopping_amplitude, atom_i=atom_i, atom_f=atom_f, orbital_i=orbital_i, orbital_f=orbital_f, hop_vector=hop_vector, spin_i=spin_i, spin_f=spin_f, add_time_reversal=add_time_reversal)
        

    def set_onsite(self, onsite, atom=None, orbital=None, spin=None, position_coordinates=None):
        if type(self._hamiltonian)==type(None):
            self._hamiltonian = np.zeros([self.n_dof, self.n_dof], dtype=COMPLEX) 
        self._hamiltonian += self._onsite_tensor(onsite=onsite, atom=atom, orbital=orbital, spin=spin, position_coordinates=position_coordinates)
        
    def set_hopping(self, hopping_amplitude, hop_vector=None, atom_i=None, atom_f=None, orbital_i=None, orbital_f=None, spin_i=None, spin_f=None, add_time_reversal=True, label=''):
        
        if type(self.kpts)==type(None):
            if type(self._hamiltonian)==type(None):
                self._hamiltonian = np.zeros([self.n_dof, self.n_dof], dtype=COMPLEX) 
            self._hamiltonian += self._hopping_tensor(hopping_amplitude=hopping_amplitude, k=None, atom_i=atom_i, atom_f=atom_f, orbital_i=orbital_i, orbital_f=orbital_f, hop_vector=hop_vector, spin_i=spin_i, spin_f=spin_f, add_time_reversal=add_time_reversal)
        else:
            if type(self._hamiltonian)==type(None):
                dimensions = [self.n_total_kpts,self.n_atoms_orbitals_spins,self.n_atoms_orbitals_spins]
                self._hamiltonian = np.zeros(dimensions)
            for i in range(self.n_total_kpts):
                k=self.flattened_kpts[:,i]
                tmp = self._hopping_tensor(hopping_amplitude=hopping_amplitude, k=k, atom_i=atom_i, atom_f=atom_f, orbital_i=orbital_i, orbital_f=orbital_f, hop_vector=hop_vector, spin_i=spin_i, spin_f=spin_f, add_time_reversal=add_time_reversal)
                self._hamiltonian[i] += tmp

        # if type(orbital_i)==type(None):
        #     orbital_i=self.atom(atom_i).orbitals
        # if type(orbital_f)==type(None):
        #     orbital_f=self.atom(atom_f).orbitals

        self.hoppings.append([atom_i,atom_f,hop_vector,add_time_reversal])
        self.hopping_labels.append(label)

    @property
    def hop_vectors(self):
        return [self.hoppings[v][2] for v in range(len(self.hoppings))]

    @property
    def hop_atoms(self):
        return [self.hoppings[v][:2] for v in range(len(self.hoppings))]

    def add_impurities(self, impurity_potential, position_coordinates, atom=None, orbital=None, spin=None, label=''):
        
        self.impurities.append([impurity_potential,position_coordinates,atom,orbital,spin])
        self.impurity_labels.append(label)

        self.set_onsite(onsite=impurity_potential, atom=atom, orbital=orbital, spin=spin, position_coordinates=position_coordinates)
        

    def set_magnetic_impurities(self, M, cell_coordinate, atom=None, orbitals=None):
        spin_matrix=np.einsum('i,ijk->jk',M,Pauli_vec)
        self.set_impurities(impurity_ampltiude=spin_matrix, cell_coordinates=cell_coordinates, orbitals=orbitals)

    ############################################################
    ################## Fourier transform #######################
    ############################################################
    # def fourier_transform(self, axes='all'):
    #     dimensions=self._pieces
    #     ft_axes=np.zeros(self.n_dimensions)
    #     if axes=='all':
    #         axes=np.arange(self.n_dimensions)
    #     elif type(axes)==int:
    #         axes=[axes]
    #     for i in axes:
    #         # check if k_space (true/false)
    #         if self.k_axes[i]:
    #             # k space => inv ft
    #             ft_axes[i]=-1
    #         else:
    #             # real space => ft
    #             ft_axes[i]=1
    #     ft=np.einsum('i,ij,ik->jk',ft_axes,xpts(dimensions),kpts(dimensions),optimize=True)
    #     ft=np.exp(1.0j*ft)
    #     ft=ft/np.sqrt(self.n_cells)
    #     ft=np.kron(np.eye(self.n_atoms_orbitals_spins),ft)

    #     if self._model=='bdg':
    #         ft = np.block([[ft,np.conj(ft)],[np.conj(ft),ft]])

    #     return ft

    # def ft_hamiltonian(self, axes='all'):
    #     ft = self.fourier_transform(axes=axes)
    #     self._hamiltonian=np.einsum('kx,xy,qy->kq',ft,self._hamiltonian,np.conjugate(ft),optimize=True)

    # def ft_wavefunction(self, axes='all'):
    #     ft = self.fourier_transform(axes=axes)
    #     self.eigenvectors = np.einsum('kx,xy->ky',ft,self.eigenvectors,optimize=True)

    ############################################################
    ################### diagonalisation ########################
    ############################################################

    def solve(self):
        t = time.time()
        if type(self.kpts)==type(None):
            self.eigenvalues,self.eigenvectors = la.eigh(self._hamiltonian, overwrite_a=True)
        elif not self.bulk_calculation:
            self._hamiltonian = scipy.linalg.block_diag(*self._hamiltonian)
            self.eigenvalues, self.eigenvectors = la.eigh(self._hamiltonian, overwrite_a=True)
        else:
            dim = self.n_atoms_orbitals_spins # dimension of (BdG) SO-matrix
            self.eigenvalues, self.eigenvectors = np.zeros([self.n_total_kpts,dim]), np.zeros([self.n_total_kpts,dim,dim])
            for dim in range(self.n_dimensions):
                for i in range(self.n_total_kpts):
                    self.eigenvalues[i], self.eigenvectors[i] = la.eigh(self._hamiltonian[i], overwrite_a=True)

        # For k lowest eigenvalues
        # from scipy.sparse.linalg import eigsh
        # w,v = eigsh(ham, k=50, which='SM', return_eigenvectors=True)
        # if self.k_space:
        #     v=self.inv_fourier_transform_psi(v)
        #     self.k_space=False

        self.exec_time = time.time() - t

        self._reshape_eigenvectors()
        return self.eigenvalues,self.eigenvectors

    def _reshape_eigenvectors(self):
        '''7-dimensional data set: [normal/anom, x, y, z, orbital, spin, eigenvalues]'''
        n=len(self.eigenvalues)

        if np.shape(self.eigenvectors)!=(n,n):
            return self.eigenvectors

        if self.bulk_calculation:
            self.eigenvectors = np.moveaxis(self.eigenvectors,0,-1)
            self.eigenvectors = np.reshape(self.eigenvectors,self.n_kpts+[self.n_atoms_orbitals,self.n_spins,2*n], 'F')
        else:

            dim=np.append(self._extended_dimensions,[2*self.n_dof])

            self.eigenvectors = np.array([np.reshape(self.eigenvectors[:self.n_dof], dim, 'F'),np.reshape(self.eigenvectors[self.n_dof:], dim, 'F')])
  
    @property
    def axes(self):
        if np.all(self._k_axes):
            return "Dimensions are k-space."
        elif np.all(self._r_axes):
            return "Dimensions are real-space."
        else:
            return f"k-space axes are {self._momentum_axes} and real-space axes are {self._real_axes}."
        
    @property
    def _momentum_axes(self):
        axes=np.arange(self.n_dimensions)
        axes=axes[self._k_axes]
        return axes

    @property
    def _r_axes(self):
        return np.logical_not(self._k_axes)

    @property
    def _real_axes(self):
        axes=np.arange(self.n_dimensions)
        axes=axes[self._r_axes]
        return axes

    def set_k_axes(self, axes):
        if type(axes)==type(None):
            axes=[]
        for axis in axes:
            if self._k_axes[axis]:
                pass
            else:
                norm = 1/np.sqrt(self._pieces[axis])
                self.eigenvectors = norm*np.fft.fft(self.eigenvectors,axis=axis+1)
                self._k_axes[axis] = not self._k_axes[axis]

    def set_real_axes(self, axes):
        if type(axes)==type(None):
            axes=[]
        for axis in axes:
            if not self._k_axes[axis]:
                pass
            else:
                norm = 1/np.sqrt(self._pieces[axis])
                self.eigenvectors = norm*np.fft.ifft(self.eigenvectors,axis=axis+1)
                self._k_axes[axis] = not self._k_axes[axis]

###################################################################
######################## Plotting ################################
###################################################################

    def lattice(self,atom,dimensions=None):
        if type(dimensions)==type(None):
            dimensions=self._pieces
        centre=np.array(np.array(dimensions)/2,dtype=int)
        lattice=np.indices(np.array(dimensions)).T -centre
        lattice=np.dot(lattice,self.lattice_vectors)
        lattice=lattice+self.atom(atom).position
        return lattice.T

    def plot_unit_cell(self, fig, ax, atoms='all', s=1, x_cells=3, y_cells=3, include_basis_vec=True, include_hoppings=True, include_atomic_labels=True):

        if atoms=='all':
            atoms=np.arange(self.n_atoms)
        
        vec=np.array(self.lattice_vectors)
        b0=vec[0]
        b1=vec[1]
        x=np.max(np.abs(b0))
        y=np.max(np.abs(b1))
        centre=b0+b1
        
        # unit cell:
        # rectangle = plt.Rectangle(-centre/2+self.com, centre[0], centre[1], fc='lightblue',ec="black")
        # plt.gca().add_patch(rectangle)
        
        n = len(atoms)
        colors = cm.rainbow(np.linspace(0, 1, n))

        # lattice:
        for i,atom in enumerate(atoms):
            x,y=self.lattice(atom,dimensions=[x_cells,y_cells])
            ax.scatter(x,y,s=s,color=colors[i],alpha=0.3)

        for i,atom in enumerate(atoms):
            x,y=self.lattice(atom,dimensions=[1,1])
            ax.scatter(x,y,s=s,color=colors[i])
            
            # atomic labels:
            if include_atomic_labels:
                pos=self.atom(atom).position
                pos=pos-np.array([-0.05,0.07])
                label=self.atom(atom).name
                ax.annotate(text=label, xytext=pos,xy=pos)

        # basis vectors:
        shift=0.00
        if include_basis_vec:
            for i in range(self.n_dimensions): # Two basis vectors are plotted in the plane 
               ax.annotate(text='', xytext=self.lattice_vectors[i]+np.array([-shift,-shift])-centre,xy=[-shift,-shift]-centre,arrowprops=dict(arrowstyle='<-', lw=2))
               ax.annotate(text=f'$b_{i}$', xytext=0.5*self.lattice_vectors[i]+[0.03,0.03]-centre,xy=[0,0]-centre)

        # hoppings
        if include_hoppings:
            for i,hopping in enumerate(self.hoppings):
                A=self.atom(hopping[0]).position
                B=self.atom(hopping[1]).position
                B=B+np.dot(hopping[2],self.lattice_vectors)
                if hopping[3]:
                    ax.annotate(text='', xytext=B,xy=A,arrowprops={'arrowstyle': '<->', 'ls': 'dashed'})
                else:
                    ax.annotate(text='', xytext=B,xy=A,arrowprops={'arrowstyle': '<-', 'ls': 'dashed'})
                ax.annotate(text=self.hopping_labels[i], xytext=0.5*(B+A)+[0.05,0.05],xy=[0,0])
        
        # axes labels
        ax.set_xlabel(r'$x/|b_0|$')
        ax.set_ylabel(r'$y/|b_1|$')
        
        return fig,ax

    def plot_lattice(self, fig, ax, energy=None, atoms='all', plot_ldos=False, plot_magnetism=False, s=1):

        if atoms=='all':
            atoms=np.arange(self.n_atoms)

        cmax=[0]
        cmin=[0]
        
        if plot_magnetism:
            cmap=cm.seismic
            label=r'$\langle\hat{\text{M}}\rangle$'
            cmax=1
            cmin=-1

        elif plot_ldos:
            cmap=cm.YlOrRd
            label=rf'$\hat{{\mathcal{{G}}}}(\omega={energy},\mathbf{{r}})$'
            cmax=1
            cmin=0
            
        if plot_magnetism:
            for atom in atoms:
                x,y=self.lattice(atom)
                ldos_up = self.local_density_of_states(energy=energy, atom=atom, orbital=None, spin='up')
                ldos_dn = self.local_density_of_states(energy=energy, atom=atom, orbital=None, spin='dn')
                ldos=ldos_up-ldos_dn
                ldos=np.fft.fftshift(ldos,0)
                sc = ax.scatter(x,y,s=s, c=ldos, cmap=cmap, vmin=cmin, vmax=cmax)

        elif plot_ldos:
            for atom in atoms:
                x,y=self.lattice(atom)
                ldos = self.local_density_of_states(energy=energy, atom=atom, orbital=None, spin=None)
                ldos=np.fft.fftshift(ldos,0)
                sc = ax.scatter(x,y,s=s, c=ldos, cmap=cmap, vmin=cmin, vmax=cmax)

        else:
            for atom in atoms:
                x,y=self.lattice(atom)
                sc = ax.scatter(x,y,s=s)

        x=np.max(np.abs(self.lattice_vectors[:,0]))
        y=np.max(np.abs(self.lattice_vectors[:,1]))
        ax.set_xlim(-0.5*x,x*(self._pieces[0]+0.5))
        ax.set_ylim(-0.5*y,y*(self._pieces[1]+0.5))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='2%', pad=0.05)
        cbar=fig.colorbar(sc,cax=cax)

        ax.annotate(label,
                    xy=(1, 1), xycoords='axes fraction',
                    xytext=(-2, -2), textcoords='offset pixels',
                    horizontalalignment='right',
                    verticalalignment='top')

        # hoppings
        for x in range(self._pieces[0]):
            for y in range(self._pieces[1]):
                for i,hopping in enumerate(self.hoppings):
                    A=self.atom(hopping[0]).position
                    B=self.atom(hopping[1]).position
                    B=B+np.dot(hopping[2],self.lattice_vectors) 
                    centre=[x,y]
                    centre=np.dot(centre,self.lattice_vectors)
                    A=A+centre
                    B=B+centre
                    ax.annotate(text='', xytext=B,xy=A,arrowprops={'arrowstyle': '-'})
        
        return fig,ax

    def plot_band_structure_path(self, fig, ax, dos):
        xmin,xmax=np.min(self.kpts),np.max(self.kpts)
        ymin,ymax=self.emin,self.emax
        im=ax.imshow(dos.T,extent=[xmin,xmax,ymin,ymax], origin='lower', cmap='Blues')
        ax.set_title(r'Tightbinding model $\overline{{\hat{{\mathcal{{G}}}}(\omega-\mu, \mathbf{{k}})}}$')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar=fig.colorbar(im,cax=cax)
        return fig, ax

    def plot_band_structure_2D(self, fig, ax, dos):
        xmin,xmax=np.min(self.kpts[0]),np.max(self.kpts[0])
        ymin,ymax=np.min(self.kpts[1]),np.max(self.kpts[1])
        im=ax.imshow(dos.T,extent=[xmin,xmax,ymin,ymax], origin='lower', cmap='Blues')
        ax.set_title(r'Tightbinding model $\overline{{\hat{{\mathcal{{G}}}}(\omega-\mu, \mathbf{{k}})}}$')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar=fig.colorbar(im,cax=cax)
        return fig, ax

    def plot_band_structure_3D(self, fig, ax, dos):

        import plotly.graph_objects as go

        X,Y,Z=self.kpts
        values=dos
        minv=np.min(values)
        maxv=np.max(values)

        fig = go.Figure(data=go.Volume(
            x=X.flatten(),
            y=Y.flatten(),
            z=Z.flatten(),
            value=values.flatten(),
           isomin=maxv/4,
           isomax=maxv,
           opacity=1, # needs to be small to see through all surfaces
           opacityscale='extremes',
           surface_count=10, # needs to be a large number for good volume rendering
           caps= dict(x_show=False, y_show=False, z_show=False), # no caps
           ))

        x=np.pi
        y=x
        z=x

        fig.update_layout(
           scene = dict(
               xaxis = dict(nticks=3, tickvals=[-x,0,x],),
               yaxis = dict(nticks=3, tickvals=[-y,0,y],),
               zaxis = dict(nticks=3, tickvals=[-z,0,z],),))

        fig.show()
        return fig, ax

    def band_structure(self, fig, ax, atom=None, oribital=None, spins=None):

       self.cmap = ListedColormap(cm.afmhot(np.linspace(0, 0.75, 256)))
       self.interpolation = 'none'

       self.xlabel=['$x/a_x$','$y/a_y$']
       self.ylabel='$\omega$'
        
       dos = self.data.density_of_states
       # dos = np.sum(dos,axis=self._extended_dimensions)

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

       return fig, ax

       def plot_fields(self, xaxis, xlabel, field, label, twin_field=None, twin_label=None, second_twin_field=None, second_twin_label=None):
        
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

        # for lattice_vector in self.lattice_vectors:


#        ldos = self.data.local_density_of_states(energy, orbitals)
#        ldos = ldos[layer]

#        # Colors and labels:
#        self.color_orb = plt.cm.rainbow(np.linspace(0,1,self.n_orbitals))
#        self.label_basis = list(string.ascii_uppercase)
#        self.label_orb = list(string .ascii_lowercase)

#        # Show hopping:
#        if show_hopping:
#            for i, link in enumerate(self.hoppings):
#                self._plt_hopping(*link, self.label_hoppings[i], show_only_centre)

#        # Find orbitals in layer to be plotted:
#        for o, orb in enumerate(self.orbitals):
#            tmp=layer+(o,slice(2))
#            print(np.shape(self.coordinates_to_orbital_position))
#            temp = self.coordinates_to_orbital_position[tmp]
#            print(np.shape(temp))
#            pos = temp.reshape([self.dimensions[0]*self.dimensions[1], 2])
#            size = ldos[...,o].reshape([self.dimensions[0]*self.dimensions[1]])
#            xx, yy = zip(*pos)
#            s=size*40/(np.mean(size))
#            self.ax.scatter(xx, yy, s=s, color=self.color_orb[o])

#        # Impurities
#        if show_impurities:
#            for impurities in self.impurities:
#                [impurity_loc, impurity_spin, impurity_orb] = impurities
#                ll=[]
#                for loc in impurity_loc:
#                    if len(loc)==2:
#                        loc.append(0)
#                    ll.append(loc)
#                ll = np.array(ll)
#                oo=[]
#                for orb in impurity_orb:
#                    oo.append(orb)
#                n_imp = np.shape(ll)[0]
#                for o in oo:
#                    if ll!=[]:
#                        indices = tuple([list(ll[:,i]) for i in range(self.n_dim)]+[o,slice(None,2,1)])
#                        pos = self.coordinates_to_orbital_position[indices]
#                        xx, yy = zip(*pos)
#                        s=10/np.mean(self.dimensions[0:1])
#                        self.ax.scatter(xx, yy, s=6*s, color='k', marker='o')
#                        self.ax.scatter(xx, yy, s=s, color=self.color_orb[o], marker='*')
        
#        # Cell borders:
#        if show_cell_borders==True:
#            width = self.basis[0][0]
#            # width = 2*basis[0][0] # 2 for hexagonal self
#            height = self.basis[1][1]

#            tmp=layer+(0,slice(2))
#            temp = self.coordinates_to_orbital_position[tmp]
#            for a_x, a_y in zip(*(temp.reshape([self.dimensions[0]*self.dimensions[1],2])+self.com[:2]).T):
#                self.ax.add_patch(Rectangle(
#                    xy=(a_x-width/2, a_y-height/2) ,width=width, height=height,
#                    linewidth=1, color='blue', fill=False))
        
#        # Annotate orbitals:
#        for o, orb in enumerate(self.orbitals):
#            if annotate_orbs==True:
#                self.ax.annotate(self.label_orb[o], self.orbitals[o][:2]+np.array([0.02,0.02]), horizontalalignment='right', verticalalignment='bottom', color=self.color_orb[o])

#        # Draw basis vectors:
#        if show_basis_vectors==True:
#            for i in range(2): # Two basis vectors are plotted in the plane
#                self.ax.annotate(text='', xytext=self.basis[i][:2],xy=self.orbitals[0][:2],  arrowprops=dict(arrowstyle='<-', lw=2))
#                self.ax.annotate(text=self.label_basis[i], xytext=0.5*self.basis[i][:2]+[0.05,0.05],xy=self.orbitals[0][:2])

#        return self.fig, self.ax

#    def _plt_hopping(self, orb_i: int, orb_f: int, cell_hop: tuple, label='', PBC=True, show_only_centre=True):
#        cell_hop = np.array(cell_hop)
#        cells = self.n_cells
#        if (cell_hop==np.zeros([self.n_dim])).all() and PBC==True:
#            if show_only_centre:
#                index=[0]
#            else:
#                index=np.arange(cells)
#            for i in index:
#                index_i = self.coord[i][:self.n_dim]
#                index_f = np.mod(index_i + cell_hop, self.dimensions)
#                index_i = list(index_i)+[orb_i]
#                orb_0 = self.coordinates_to_orbital_position[tuple(index_i)] 
#                index_f = list(index_f)+[orb_f]
#                orb_1 = self.coordinates_to_orbital_position[tuple(index_f)]
#                orb_0=orb_0[:2]
#                orb_1=orb_1[:2]
#                self.ax.annotate(text='', xytext=orb_1, xy=orb_0,  arrowprops={'arrowstyle': '<->', 'ls': 'dashed'})
#                self.ax.annotate(text=label, xytext=0.5*(orb_0+orb_1)+(0,0.01), xy=orb_0)
#        # else:
#        for dimension in range(self.n_dim):
#            temp=np.zeros(self.n_dim)
#            temp[dimension]=1
#            if (temp==cell_hop).all():
#                if show_only_centre:
#                    index=[0]
#                else:
#                    index=np.arange(cells)
#                for i in index: 
#                    index_i = self.coord[i][:self.n_dim]
#                    index_f = np.mod(index_i + cell_hop, self.dimensions)
#                    coord_i = self.coord_cell[i]
#                    edge = bool(coord_i[dimension]>=self.centre[dimension])
#                    if edge and not PBC:
#                        pass
#                    else:
#                        index_i = list(index_i)+[orb_i,ALL]
#                        orb_0 = self.coordinates_to_orbital_position[tuple(index_i)] 
#                        index_f = list(index_f)+[orb_f,ALL]
#                        orb_1 = self.coordinates_to_orbital_position[tuple(index_f)]
#                        orb_0=orb_0[:2]
#                        orb_1=orb_1[:2]
#                        self.ax.annotate(text='', xytext=orb_1, xy=orb_0,  arrowprops={'arrowstyle': '<->', 'ls': 'dashed'})
#                        self.ax.annotate(text=label, xytext=0.5*(orb_0+orb_1)+(0,0.01), xy=orb_0)




###################################################################
##################### BdG model ##########################
###################################################################

class BogoliubovdeGennes(Tightbinding):
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

    def __init__(self,lattice_vectors,name):

        super().__init__(lattice_vectors,name)

        self._model='bdg'

        self.iterations=0

        self._hartree=None
        self._fock=None
        self._gorkov=None

        self._hartree_indices = []
        self._fock_indices = []
        self._gorkov_indices = []

        self._hubbard_u = None
        
        self._hartree_print_indices = []
        self._fock_print_indices = []
        self._gorkov_print_indices = []

        self._hartree_iterations=[]
        self._fock_iterations=[]
        self._gorkov_iterations=[]
        
        self.V=[]
        self.V_mf=[]
        self.Eg=[]
        self.f_mf=[]
        self.free_energy=[]
        
        self.print_V=False
        self.print_V_mf=False
        self.print_Eg=False
        self.print_free_energy=False
    
    def set_friction(friction=0.7):
        self.friction = 0.7

    def set_max_iterations(self,max_iterations=100):
        self.max_iterations = 1000

    def set_absolute_convergence_factor(self,absolute_convergence_factor=0.00001):
        self.absolute_convergence_factor = absolute_convergence_factor

    def _set_mean_field_hamiltonian(self):

        # initialise tight binding hamiltonian if doesn't exist:
        try: 
            self._tb_ham
        except:
            self._tb_ham = self._hamiltonian

        n_dof = self.n_dof

        try:
            self._hubbard_indices
        except:
            self._set_hubbard_indices()

            if type(self._hartree)==type(None):
                self._hartree = np.zeros(self.n_dof, dtype=COMPLEX) 
            if type(self._fock)==type(None):
                self._fock = np.zeros([self.n_dof,self.n_dof], dtype=COMPLEX) 
            if type(self._gorkov)==type(None):
                self._gorkov = np.zeros([self.n_dof,self.n_dof], dtype=COMPLEX) 
            self._fock=self._fock[self._hubbard_indices]
            self._gorkov=self._gorkov[self._hubbard_indices]

        hartree=self._hartree
        fock=self._fock
        gorkov=self._gorkov

        hamiltonian = np.zeros([2*n_dof,2*n_dof],dtype=COMPLEX)
        hamiltonian[:n_dof,:n_dof]=self._tb_ham-np.diag(hartree)
        hamiltonian[n_dof:,n_dof:]=-(np.conj(self._tb_ham)-np.conj(np.diag(hartree)))
        hamiltonian[self._hubbard_indices[0],self._anomalous_indices[1]]=-fock
        hamiltonian[self._anomalous_indices[0],self._hubbard_indices[1]]=-(-np.conj(fock))
        hamiltonian[self._hubbard_indices[0],self._anomalous_indices[1]]=-np.conj(-gorkov)
        hamiltonian[self._anomalous_indices[0],self._hubbard_indices[1]]=-gorkov
        self._hamiltonian = hamiltonian

    def reset_hartree(self):
        self._hartree=np.zeros([self.n_dof], dtype=COMPLEX)
        self._hartree_iterations=[]

    def reset_fock(self):
        self._fock=np.zeros([self.n_dof,self.n_dof], dtype=COMPLEX)
        self._fock_iterations=[]

    def reset_gorkov(self):
        self._gorkov=np.zeros([self.n_dof,self.n_dof], dtype=COMPLEX)
        self._gorkov_iterations=[]

        self._hartree=np.zeros([self.n_dof], dtype=COMPLEX)
        self._fock=np.zeros([self.n_dof,self.n_dof], dtype=COMPLEX)
        self._gorkov=np.zeros([self.n_dof,self.n_dof], dtype=COMPLEX)

    def set_hartree(self,onsite,atom=None,orbital=None,spin=None,position_coordinates=None):
        if type(self._hartree)==type(None):
            self._hartree = np.zeros(self.n_dof, dtype=COMPLEX) 
        self._hartree += np.diag(self._onsite_tensor(onsite=onsite, atom=atom, orbital=orbital, spin=spin, position_coordinates=position_coordinates))

    def set_fock(self, hopping_amplitude, k=None, atom_i=None, atom_f=None, orbital_i=None, orbital_f=None, spin_i=None, spin_f=None, hop_vector=None, add_time_reversal=True):
        if type(self._fock)==type(None):
            self._fock = np.zeros([self.n_dof,self.n_dof], dtype=COMPLEX) 
        self._fock += self._hopping_tensor(hopping_amplitude=hopping_amplitude, k=k, atom_i=atom_i, atom_f=atom_f, orbital_i=orbital_i, orbital_f=orbital_f, hop_vector=hop_vector, spin_i=spin_i, spin_f=spin_f, add_time_reversal=add_time_reversal)

    def set_gorkov(self, hopping_amplitude, k=None, atom_i=None, atom_f=None, orbital_i=None, orbital_f=None, spin_i=None, spin_f=None, hop_vector=None, add_time_reversal=True):
        if type(self._gorkov)==type(None):
            self._gorkov = np.zeros([self.n_dof,self.n_dof], dtype=COMPLEX) 
        self._gorkov += self._hopping_tensor(hopping_amplitude=hopping_amplitude, k=k, atom_i=atom_i, atom_f=atom_f, orbital_i=orbital_i, orbital_f=orbital_f, hop_vector=hop_vector, spin_i=spin_i, spin_f=spin_f, add_time_reversal=add_time_reversal)

    def set_external_hartree(self,onsite,atom=None,orbital=None,spin=None,position_coordinates=None):
        self.external_hartree += np.diag(self._onsite_tensor(onsite=onsite, atom=atom, orbital=orbital, spin=spin, position_coordinates=position_coordinates))

    def set_external_fock(self, hopping_amplitude, k=None, atom_i=None, atom_f=None, orbital_i=None, orbital_f=None, spin_i=None, spin_f=None, hop_vector=None, add_time_reversal=True):
        self.external_fock += self._hopping_tensor(hopping_amplitude=hopping_amplitude, k=k, atom_i=atom_i, atom_f=atom_f, orbital_i=orbital_i, orbital_f=orbital_f, hop_vector=hop_vector, spin_i=spin_i, spin_f=spin_f, add_time_reversal=add_time_reversal)

    def set_external_gorkov(self, hopping_amplitude, k=None, atom_i=None, atom_f=None, orbital_i=None, orbital_f=None, spin_i=None, spin_f=None, hop_vector=None, add_time_reversal=True):
        self.external_gorkov += self._hopping_tensor(hopping_amplitude=hopping_amplitude, k=k, atom_i=atom_i, atom_f=atom_f, orbital_i=orbital_i, orbital_f=orbital_f, hop_vector=hop_vector, spin_i=spin_i, spin_f=spin_f, add_time_reversal=add_time_reversal)

    def set_hubbard_u(self, hopping_amplitude, k=None, atom_i=None, atom_f=None, orbital_i=None, orbital_f=None, spin_i=None, spin_f=None, hop_vector=None, add_time_reversal=True):
        if type(self._hubbard_u)==type(None):
            self._hubbard_u = np.zeros([self.n_dof, self.n_dof], dtype=COMPLEX) 
        self._hubbard_u += self._hopping_tensor(hopping_amplitude=hopping_amplitude, k=k, atom_i=atom_i, atom_f=atom_f, orbital_i=orbital_i, orbital_f=orbital_f, hop_vector=hop_vector, spin_i=spin_i, spin_f=spin_f, add_time_reversal=add_time_reversal)

    def _set_hubbard_indices(self):

        self._hubbard_indices = np.nonzero(self._hubbard_u)

        self.U_entries = self._hubbard_u[self._hubbard_indices]

        self._anomalous_indices=np.array(self._hubbard_indices)
        self._anomalous_indices+=self.n_dof
        self._anomalous_indices=tuple(self._anomalous_indices)

    def record_hartree(self, location, atom, orbital=None, spin=None, _print=False):
        if type(orbital)==type(None):
            orbital=0
        if type(spin)==type(None) and self.n_spins==1:
            spin=0
        orbital=self.index(atom,orbital)
        spin=self.spin(spin)
        self._hartree_print = _print
        tmp=np.append(np.append(location,orbital),spin)
        index=np.ravel(coordinates_to_indices(tmp, self._extended_dimensions))
        self._hartree_indices.append(index)

    def record_fock(self, location_i, location_f, atom_i, atom_f, orbital_i=None, orbital_f=None, spin_i=None, spin_f=None, _print=False):
        self._fock_print = _print
        if type(orbital_i)==type(None) and type(orbital_f)==type(None):
            orbital_i=orbital_f=0
        if type(spin_i)==type(None) and type(spin_f)==type(None) and self.n_spins==1:
            spin_i=spin_f=0
        orbital_i=self.index(atom_i,orbital_i,spin=0)
        orbital_f=self.index(atom_f,orbital_f,spin=0)
        tmp_i=np.append(np.append(location_i,orbital_i),spin_i)
        tmp_f=np.append(np.append(location_f,orbital_f),spin_f)
        index_i=coordinates_to_indices(tmp_i, self._extended_dimensions)
        index_f=coordinates_to_indices(tmp_f, self._extended_dimensions)
        temp=np.zeros([self.n_dof,self.n_dof])
        temp[index_i,index_f]=1
        temp = np.ravel(np.nonzero(temp))
        self._fock_indices.append(temp)

    def record_gorkov(self, location_i, location_f, atom_i, atom_f, orbital_i=None, orbital_f=None, spin_i=None, spin_f=None, _print=False):
        self._gorkov_print = _print
        if type(orbital_i)==type(None) and type(orbital_f)==type(None):
            orbital_i=orbital_f=0
        if type(spin_i)==type(None) and type(spin_f)==type(None) and self.n_spins==1:
            spin_i=spin_f=0
        orbital_i=self.index(atom_i,orbital_i,spin=0)
        orbital_f=self.index(atom_f,orbital_f,spin=0)
        tmp_i=np.append(np.append(location_i,orbital_i),spin_i)
        # tmp_i=np.append(np.append(location_i,spin_i),orbital_i)
        tmp_f=np.append(np.append(location_f,orbital_f),spin_f)
        # tmp_f=np.append(np.append(location_f,spin_f),orbital_f)
        index_i=coordinates_to_indices(tmp_i, self._extended_dimensions)
        index_f=coordinates_to_indices(tmp_f, self._extended_dimensions)
        temp=np.zeros([self.n_dof,self.n_dof])
        temp[index_i,index_f]=1
        temp = np.ravel(np.nonzero(temp))
        self._gorkov_indices.append(temp)

    def set_friction(self,friction):
        self.friction = friction

    def set_max_iterations(self,max_iterations):
        self.max_iterations = max_iterations

    def set_absolute_convergence_factor(self,absolute_convergence_factor):
        self.absolute_convergence_factor = absolute_convergence_factor

    #########################################
    ########## Statistical Mechanics ########
    #########################################
    
    @property
    def anomalous_density_matrix(self):
        return np.einsum('in,in->in',np.roll(self.eigenvectors[:self.n_dof],self.n_dof,axis=0),np.conj(self.eigenvectors[:self.n_dof]),optimize=True)

    @property
    def thermal_density_matrix(self):
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
        w,v=self.eigenvalues,self.eigenvectors
        T=self.temperature
        tmp0=v[self._hubbard_indices[0]]
        tmp1=np.conj(v)[self._hubbard_indices[1]]
        tmp1a=v[self._anomalous_indices[1]]
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

    def _set_fields(self, trace_density, density, anomalous_density):
        """Returns Hartree, Fock, Gorkov"""
        hartree = +np.einsum('ij,j->i',self._hubbard_u,trace_density,optimize=True)
        fock    = -np.multiply(self.U_entries,density)
        gorkov  = np.multiply(self.U_entries,anomalous_density)
        return np.real_if_close(hartree), np.real_if_close(fock), np.real_if_close(gorkov)

    def calculate_free_energy(self, trace_density, density, anomalous_density):

        self.Eg.append(np.real(np.sum(self.eigenvalues[:self.n_dof] + np.diagonal(self._hamiltonian)[:self.n_dof])/self.n_dof))

        if self.temperature==0:
            self.f_mf.append(self.Eg[-1])
        else:
            self.f_mf.append(np.real(self.Eg[-1]-2*self.temperature*np.sum(np.log(1+np.exp(-self.eigenvalues[self.n_dof:]/self.temperature)))/self.n_dof))
        
        h=np.sum(self._hartree*trace_density)
        h+=h
        f=np.sum(self._fock*density)
        f+=np.conj(f)
        g=np.sum(self._gorkov*anomalous_density)
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
            temp_hartree.append(complex(self._hartree[index]))
        
        for index in self._fock_indices:
            for i,hubbard_index in enumerate(hubbard_indices):
                if np.array_equal(index,hubbard_index):
                    temp_fock.append(self._fock[i])

        for index in self._gorkov_indices:
            for i,hubbard_index in enumerate(hubbard_indices):
                if np.array_equal(index,hubbard_index):
                    temp_gorkov.append(self._gorkov[i])

        self._hartree_iterations.append(temp_hartree)
        self._fock_iterations.append(temp_fock)
        self._gorkov_iterations.append(temp_gorkov)

        if np.sum(self._hartree_print_indices)>0:
            print('\nHartree field:')
            temp=np.array(self._hartree_iterations[-1])
            print(temp[self._hartree_print_indices])
        if np.sum(self._fock_print_indices)>0:
            print('\nFock field:')
            temp=np.array(self._fock_iterations[-1])
            print(temp[self._fock_print_indices])
        if np.sum(self._gorkov_print_indices)>0:
            print('\nGorkov field:')
            temp=np.array(self._gorkov_iterations[-1])
            print(temp[self._gorkov_print_indices])

    def solve(self,reshape=True):
        """Decorates self.solve() with self._set_mean_field_hamiltonian()"""

        self._set_mean_field_hamiltonian()

        t = time.time()
        self.eigenvalues,self.eigenvectors = la.eigh(self._hamiltonian, overwrite_a=True)
        self.exec_time = time.time() - t
        if reshape:
            self._reshape_eigenvectors()
        return self.eigenvalues,self.eigenvectors

    def __iter__(self):

        return self

    def __next__(self):
        
        ################## Record fields ###################
        
        self.record_fields()

        ############ Solve mean-field Hamiltonian #############

        self.solve(reshape=False)

        trace_density, density, anomalous_density = self.thermal_density_matrix

        hartree, fock, gorkov = self._set_fields(trace_density, density, anomalous_density)

        ##################### Free energy ######################
        
        self.calculate_free_energy(trace_density, density, anomalous_density)

        ####################### Friction ########################

        hartree = (1-self.friction)*hartree + self.friction*self._hartree
        fock = (1-self.friction)*fock + self.friction*self._fock
        gorkov = (1-self.friction)*gorkov + self.friction*self._gorkov
        
        ###################### Update fields #####################

        self._hartree, self._fock, self._gorkov = hartree, fock, gorkov 

        return self.eigenvalues,self.eigenvectors

    def self_consistent_calculation(self,friction=0.7,max_iterations=100,absolute_convergence_factor=0.00001):
        """If dos=True, the density of states are calculated once the self consistent
        loop has converged."""

        self.set_friction(friction=friction)
        self.set_max_iterations(max_iterations=max_iterations)
        self.set_absolute_convergence_factor(absolute_convergence_factor=absolute_convergence_factor)

        t = time.time()

        self._set_mean_field_hamiltonian()

        iteration = iter(self)
        
        for i in tqdm.tqdm(range(self.max_iterations)):
        # for i in range(self.max_iterations):
            
            self.iterations+=1

            hartree_old = self._hartree
            fock_old = self._fock
            gorkov_old = self._gorkov

            next(iteration)
            
            eps = self.absolute_convergence_factor
            if (np.allclose(hartree_old,self._hartree,atol=eps) & np.allclose(fock_old,self._fock,atol=eps) & np.allclose(gorkov_old,self._gorkov,atol=eps)):
                self.converged=True
                break
            if i+1 == self.max_iterations:
                self.converged=False
                print('Did not converge within max_iterations!')
                break

        self.exec_time = time.time() - t

        self._hartree_iterations = np.real_if_close(self._hartree_iterations).T
        self._fock_iterations = np.real_if_close(self._fock_iterations).T
        self._gorkov_iterations = np.real_if_close(self._gorkov_iterations).T
        self._reshape_eigenvectors()
    
    @property
    def hartree_array(self):
        temp=np.reshape(self._hartree,self._extended_dimensions,'F')
        return np.real_if_close(temp)

    @property
    def fock_array(self):
        extended_dimensions=np.append(self._extended_dimensions,self._extended_dimensions)
        temp=np.zeros([self.n_dof,self.n_dof],dtype=COMPLEX)
        temp[self._hubbard_indices]=self._fock
        temp=np.reshape(temp,extended_dimensions,'F')
        return np.real_if_close(temp)
    
    @property
    def gorkov_array(self):
        extended_dimensions=np.append(self._extended_dimensions,self._extended_dimensions)
        temp=np.zeros([self.n_dof,self.n_dof],dtype=COMPLEX)
        temp[self._hubbard_indices]=self._gorkov
        temp=np.reshape(temp,extended_dimensions,'F')
        return np.real_if_close(temp)

    def hartree(self,atom=None,orbital=None,spin=None):
        
        indices=[]
        if type(atom)!=type(None):
            _atoms=[self._atom_index(atom)]
        else:
            _atoms=self._atom_indices
        for atom in _atoms:
            n_orbitals=self.n_orbitals(atom)
            if type(orbital)==type(None):
                for orb in np.arange(self.n_orbitals(atom)):
                    indices.append(self.index(atom,orb))
            else:
                indices.append(self.index(atom,orbital))
        
        hartree=np.copy(self.hartree_array)
        if type(spin)==type(None):
            hartree=np.sum(hartree,-1)
        else:
            hartree=hartree[...,spin]
        hartree=np.sum(hartree[...,indices],-1)

        return hartree

    def _field_indexing(self, field_array, atom_i=None, atom_f=None, orbital_i=None, orbital_f=None, hop_vector=[0,0], spin_i=None, spin_f=None):

        if type(hop_vector)==type(None):
            hop_vector=[0,0]

        if type(spin_i)==type(None) and type(spin_f)==type(None):
            spins_i=np.arange(self.n_spins)
            spins_f=np.arange(self.n_spins)
        elif type(spin_i)==str or type(spin_i)==int and type(spin_i)==str or type(spin_i)==int:
            spins_i=[self.spin(spin_i)]
            spins_f=[self.spin(spin_f)]
        else:
            raise ValueError('Both spin_i and spin_f must be None or both must be type either int or str!')

        if type(atom_i)==type(None) and type(atom_f)==type(None):
            atoms_i=np.arange(self.n_atoms)
            atoms_f=np.arange(self.n_atoms)
        elif type(atom_i)==str or type(atom_i)==int and type(atom_i)==str or type(atom_i)==int:
            atoms_i=[self.atom(atom_i)]
            atoms_f=[self.atom(atom_f)]
        else:
            raise ValueError('Both atom_i and atom_f must be None or both must be type either int or str!')


        if type(orbital_i)==type(None) and type(orbital_f)==type(None):
            if type(atom_i)==type(None) and type(atom_f)==type(None):
                indices_i=np.arange(self.n_atoms_orbitals)
                indices_f=indices_i
            else:
                
                orbitals_i=np.arange(self.atom(atom_i).n_orbitals)
                orbitals_f=np.arange(self.atom(atom_f).n_orbitals)
                indices_i=[self.index(atom_i,o) for o in orbitals_i]
                indices_f=[self.index(atom_f,o) for o in orbitals_f]
        elif type(atom_i)==str or type(atom_i)==int and type(atom_i)==str or type(atom_i)==int:
                indices_i=[self.index(a,orbital_i) for a in atoms_i]
                indices_f=[self.index(a,orbital_f) for a in atoms_f]
        else:
            raise ValueError('Both orbital_i and orbital_f must be None or both must be type either int or str!')
        
        l=self.n_dimensions+2
        for axis in range(len(hop_vector)):
            field_array=np.roll(field_array,hop_vector[axis],axis=l+axis)

        for i in range(self.n_dimensions):
            field_array=np.diagonal(field_array,axis1=0,axis2=l-i)
        for i in range(self.n_dimensions):
            field_array=np.moveaxis(field_array, -1, 0)

        field_array=np.array([[field_array[...,indices_f[j],spins_f[i],indices_i[j],spins_i[i]] for i in range(len(spins_f))] for j in range(len(indices_i))] )

        field_array=np.sum(field_array, 0)/len(indices_i)
        field_array=np.sum(field_array, 0)/len(spins_i)

        return field_array

    def fock(self, atom_i=None, atom_f=None, orbital_i=None, orbital_f=None, hop_vector=None, spin_i=None, spin_f=None):
        return self._field_indexing(field_array=self.fock_array, atom_i=atom_i, atom_f=atom_f, orbital_i=orbital_i, orbital_f=orbital_f, hop_vector=hop_vector, spin_i=spin_i, spin_f=spin_f)

    def gorkov(self, atom_i=None, atom_f=None, orbital_i=None, orbital_f=None, hop_vector=None, spin_i=None, spin_f=None):
        return self._field_indexing(field_array=self.gorkov_array, atom_i=atom_i, atom_f=atom_f, orbital_i=orbital_i, orbital_f=orbital_f, hop_vector=hop_vector, spin_i=spin_i, spin_f=spin_f)
###################################################################
##################### Plotting ##########################
###################################################################
        
###################################################################


#class plot_data(processed_data):

#    def __init__(self, fig, ax, data):
        
#        self.fig = fig
#        self.ax = ax
        
#        self.data=data
#        self.energy_interval=data.energy_interval
#        self.resolution=data.resolution
#        self.density_of_states=data.density_of_states
#        self.n_spins = data.n_spins
#        self.n_dim = data.n_dim
#        self.hoppings = data.hoppings
#        self.label_hoppings = data.label_hoppings
#        self.impurities = data.impurities
#        self.n_orbitals = data.n_orbitals
#        self.orbitals = data.orbitals
#        self.n_cells = data.n_cells
#        self.dimensions = data.dimensions
#        self.coord = data.coord
#        self.coord_cell = data.coord_cell
#        self.centre = data.centre
#        self.coordinates_to_orbital_position = data.coordinates_to_orbital_position
#        self.basis = data.basis
#        self.com = data.com

#        self.hartree_iterations=data.hartree_iterations
#        self.fock_iterations=data.fock_iterations
#        self.gorkov_iterations=data.gorkov_iterations
#        self.V=data.V
#        self.V_mf=data.V_mf
#        self.Eg=data.Eg
#        self.free_energy=data.free_energy

#    def _imshow(self, ldos):

#        self.interpolation = 'none'

#        x=np.fft.fftshift(ldos)

#        self.im = self.ax.imshow(
#                x[::-1].T, extent=self._extent,
#                interpolation=self.interpolation,
#                vmin=self.vmin, vmax=self.vmax,
#                cmap=self.cmap,
#                origin='lower')

#    def differential_current_map(self, energy, layer=(ALL,ALL,0), orbital=0, spin=None, cartesian=False, n_pts=10, gaussian_mean=0.4):

#        x,y=self.dimensions[:2]
#        self._extent = [-x/2,x/2,-y/2,y/2]

#        ldos = self.data.local_density_of_states(energy, orbital, spin)

#        ldos = ldos[layer]
#        if cartesian:
#            ldos=self.data.ldos_cartesian(ldos, n_pts, gaussian_mean)
#            self._extent = [-x/2-0.5/n_pts,x/2+0.5/n_pts,-y/2-0.5/n_pts,y/2+0.5/n_pts]

#        eV = energy
#        epsilon = self.resolution

#        self.text= (f'DOS map'
#        '\n'
#        f'$\omega={eV:.2f}$'
#        '\n'
#        f'$\epsilon={epsilon}$'
#        )    

#        self.vmin, self.vmax = np.amin(ldos), np.amax(ldos)

#        if layer==(ALL,ALL,0):
#            self.ax.set(xlabel='$x/a_x$', ylabel='$y/a_y$')

#        if layer==(ALL,0,ALL):
#            self.ax.set(xlabel='$x/a_x$', ylabel='$z/a_z$')

#        if layer==(0,ALL,ALL):
#            self.ax.set(xlabel='$y/a_y$', ylabel='$z/a_z$')

#        self.cmap = ListedColormap(cm.afmhot(np.linspace(0, 0.75, 256)))

#        self._imshow(ldos)

#        self.fig.subplots_adjust(hspace=0.1, wspace=0.1)

#        return self.fig,self.ax

#    def probe_magnetic_z(self, energy, layer=(ALL,ALL,0), orbital=0):

#        ldos = self.data.local_density_of_states(energy, orbital)

#        ldos = ldos[...,0] - ldos[...,1]
    
#        ldos = ldos[layer]

#        eV = energy
#        epsilon = self.resolution

#        self.text= (f'\map'
#        '\n'
#        f'$\omega={eV:.2f}$'
#        '\n'
#        f'$\epsilon={epsilon}$'
#        )    

#        self.vmin, self.vmax = np.amin(ldos), np.amax(ldos)

#        x,y=np.shape(ldos)[:2]

#        if layer==(ALL,ALL,0):
#            self.ax.set(xlabel='$x/a_x$', ylabel='$y/a_y$')

#        if layer==(ALL,0,ALL):
#            self.ax.set(xlabel='$x/a_x$', ylabel='$z/a_z$')

#        if layer==(0,ALL,ALL):
#            self.ax.set(xlabel='$y/a_y$', ylabel='$z/a_z$')

#        self._extent = [-x/2,x/2,-y/2,y/2]
        
#        self.cmap = ListedColormap(cm.PiYG(np.linspace(0, 1, 256)))

#        self._imshow(ldos)

#        self.fig.subplots_adjust(hspace=0.1, wspace=0.1)

#        return self.fig,self.ax

#    def quasiparticle_interference(self, energy, layer=(ALL,ALL,0), orbital=0, spins=None, remove_central_bright_spot=True):
        
#        ldos = self.ldos=greens_function.local_density_of_states(energy, orbital, spin)
#        ldos = self.ldos[layer]

#        self.vmin, self.vmax = np.amin(ldos), np.amax(ldos)

#        f = np.fft.fft2(ldos, axes=(0,1), norm='ortho')
#        abs_f = np.abs(f)

#        if layer==(ALL,ALL,0):
#            self.ax.set(xlabel='$k_x a_x$', ylabel='$k_y a_y$')

#        if layer==(ALL,0,ALL):
#            self.ax.set(xlabel='$k_x a_x$', ylabel='$k_z a_z$')

#        if layer==(0,ALL,ALL):
#            self.ax.set(xlabel='$k_y a_y$', ylabel='$k_z a_z$')
        
#        self._extent = [-1,1,-1,1]

#        x_label_list = ['$-\pi$', '0', '$\pi$']
#        y_label_list = ['$-\pi$', '0', '$\pi$']

#        self.ax.set_xticks([-1,0,1])
#        self.ax.set_yticks([-1,0,1])

#        self.ax.set_xticklabels(x_label_list)
#        self.ax.set_yticklabels(x_label_list)


#        # max/min without central bright spot and lines:
#        if remove_central_bright_spot:
#            vmax=[]
#            vmin=[]
#            temp=np.copy(abs_f)
#            temp[0,0]=0
#            vmax.append(np.amax(temp))
#            vmin.append(np.amin(temp))
#            self.vmax = np.max(vmax)
#            self.vmin = np.max(vmin)
        
#        self.cmap = LinearSegmentedColormap.from_list("", ["black","orange","white"])

#        self._imshow(abs_f)

#        self.fig.subplots_adjust(hspace=0.1, wspace=0.1)

#        return self.fig,self.ax
        
#    def reciprocal_space_surface(self, energy, layer=(ALL,ALL,0), orbital=0, spins=None):

#        ldos = self.ldos=greens_function.local_density_of_states(energy, orbital, spin)

#        ldos = self.ldos[layer]

#        self.vmin, self.vmax = np.amin(ldos), np.amax(ldos)

#        x,y=np.shape(ldos)[:2]

#        if layer==(ALL,ALL,0):
#            self.ax.set(xlabel='$k_x a_x$', ylabel='$k_y a_y$')

#        if layer==(ALL,0,ALL):
#            self.ax.set(xlabel='$k_x a_x$', ylabel='$k_z a_z$')

#        if layer==(0,ALL,ALL):
#            self.ax.set(xlabel='$k_y a_y$', ylabel='$k_z a_z$')

#        self._extent = [-1,1,-1,1]

#        x_label_list = ['$-\pi$', '0', '$\pi$']
#        y_label_list = ['$-\pi$', '0', '$\pi$']

#        self.ax.set_xticks([-1,0,1])
#        self.ax.set_yticks([-1,0,1])

#        self.ax.set_xticklabels(x_label_list)
#        self.ax.set_yticklabels(x_label_list)
        
#        self.cmap = ListedColormap(cm.afmhot(np.linspace(0, 0.75, 256)))

#        self._imshow(ldos)

#        self.fig.subplots_adjust(hspace=0.1, wspace=0.1)

#        return self.fig,self.ax


#    def set_text_box(self):
#        text_box = AnchoredText(self.text, loc=2, pad=0.3, borderpad=0)
#        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
#        self.ax.add_artist(text_box)

#    def set_cbar(self): 
#        cax = inset_axes(self.ax,width = '5%',height = '80%',loc = 5) 
#        self.cbar = plt.colorbar(self.im, ax = self.ax,cax = cax ,format = '%.2f' ,ticks = [self.vmin,self.vmax])       
#        cax.yaxis.set_ticks_position('left')

#    def set_label(self,label):
#        text_box = AnchoredText(label, loc=3, pad=0.3, borderpad=0, frameon=False)
#        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
#        self.ax.add_artist(text_box)

#    def plot_3D(self, energy, orbital=0, spins=None):

#        import plotly.graph_objects as go

#        ldos = self.ldos=greens_function.local_density_of_states(energy, orbital, spin)

#        x=int(self.dimensions[0]/2)
#        y=int(self.dimensions[1]/2)
#        z=int(self.dimensions[2]/2)
#        X, Y, Z = np.mgrid[-x:x+1, -y:y+1, -z:z+1]
#        values = np.fft.fftshift(self.ldos)
#        #values=self.ldos
#        minv=np.min(values)
#        maxv=np.max(values)

#        fig = go.Figure(data=go.Volume(
#            x=X.flatten(),
#            y=Y.flatten(),
#            z=Z.flatten(),
#            value=values.flatten(),
#            isomin=maxv/2,
#            isomax=maxv,
#            opacity=1, # needs to be small to see through all surfaces
#            opacityscale='extremes',
#            surface_count=10, # needs to be a large number for good volume rendering
#            caps= dict(x_show=False, y_show=False, z_show=False), # no caps
#            ))

#        fig.update_layout(
#            scene = dict(
#                xaxis = dict(nticks=3, tickvals=[-x,0,x],),
#                yaxis = dict(nticks=3, tickvals=[-y,0,y],),
#                zaxis = dict(nticks=3, tickvals=[-z,0,z],),))
#            #width=700,
#            #margin=dict(r=20, l=10, b=10, t=10))

#        return fig

#    def spectrum(self, locations, orbital, spins=None):

#        spectrum = self.data.spectrum(locations, orbital, spin)

#        n = len(spectrum)
#        colors = [cm.nipy_spectral(i) for i in np.linspace(0, 1, n)]

#        for i in range(n):


#            if not spin:
#                s=0 
#                x,y=self.energy_interval,spectrum[i,s]
#                self.ax.plot(x,y, linestyle='solid', color=colors[i],label=(f'$\mathbf{{r}}={locations[i]}$, ' + r'$\uparrow$'))
#                self.ax.fill_between(x,y,0,color=colors[i], alpha=0.4)
#                s=1
#                x,y=self.energy_interval,spectrum[i,s]
#                self.ax.plot(x,y, linestyle='dotted', color=colors[i],label=(f'$\mathbf{{r}}={locations[i]}$, ' + r'$\downarrow$'))
#                self.ax.fill_between(x,y,0,color=colors[i], alpha=0.4)
#            else:
#                x,y=self.energy_interval,spectrum[i]
#                self.ax.plot(x,y, linestyle='solid', color=colors[i],label=(f'$\mathbf{{r}}={locations[i]}$'))
#                self.ax.fill_between(x,y,0,color=colors[i], alpha=0.4)


#        self.ax.set(xlabel=r'$\omega$')  
#        self.ax.set(ylabel=r'Density of states')

#        text_DOS = ('Spectrum'
#        '\n'
#        f'$\epsilon={self.resolution}$'
#        )    

#        self.fig.text(.22, .81, text_DOS,
#                 {'bbox': dict(boxstyle="square", alpha=0.5, fc="white",
#                               ec="none", pad=0.2)}, ha='center', va='center')
                               
#        legend = self.fig.legend(loc="upper left", 
#        fancybox=True, shadow=True, #prop=fontP,
#        bbox_to_anchor=(0.14,0.58,1,0.2))

#        return self.fig,self.ax

#    def lattice(self, energy, model, layer, spins, orbitals, annotate_orbs=True, show_cell_borders=True, show_basis_vectors=True, show_hopping=True, show_impurities=True, show_only_centre=True):

#        ldos = self.data.local_density_of_states(energy, orbitals)
#        ldos = ldos[layer]

#        # Colors and labels:
#        self.color_orb = plt.cm.rainbow(np.linspace(0,1,self.n_orbitals))
#        self.label_basis = list(string.ascii_uppercase)
#        self.label_orb = list(string .ascii_lowercase)

#        # Show hopping:
#        if show_hopping:
#            for i, link in enumerate(self.hoppings):
#                self._plt_hopping(*link, self.label_hoppings[i], show_only_centre)

#        # Find orbitals in layer to be plotted:
#        for o, orb in enumerate(self.orbitals):
#            tmp=layer+(o,slice(2))
#            print(np.shape(self.coordinates_to_orbital_position))
#            temp = self.coordinates_to_orbital_position[tmp]
#            print(np.shape(temp))
#            pos = temp.reshape([self.dimensions[0]*self.dimensions[1], 2])
#            size = ldos[...,o].reshape([self.dimensions[0]*self.dimensions[1]])
#            xx, yy = zip(*pos)
#            s=size*40/(np.mean(size))
#            self.ax.scatter(xx, yy, s=s, color=self.color_orb[o])

#        # Impurities
#        if show_impurities:
#            for impurities in self.impurities:
#                [impurity_loc, impurity_spin, impurity_orb] = impurities
#                ll=[]
#                for loc in impurity_loc:
#                    if len(loc)==2:
#                        loc.append(0)
#                    ll.append(loc)
#                ll = np.array(ll)
#                oo=[]
#                for orb in impurity_orb:
#                    oo.append(orb)
#                n_imp = np.shape(ll)[0]
#                for o in oo:
#                    if ll!=[]:
#                        indices = tuple([list(ll[:,i]) for i in range(self.n_dim)]+[o,slice(None,2,1)])
#                        pos = self.coordinates_to_orbital_position[indices]
#                        xx, yy = zip(*pos)
#                        s=10/np.mean(self.dimensions[0:1])
#                        self.ax.scatter(xx, yy, s=6*s, color='k', marker='o')
#                        self.ax.scatter(xx, yy, s=s, color=self.color_orb[o], marker='*')
        
#        # Cell borders:
#        if show_cell_borders==True:
#            width = self.basis[0][0]
#            # width = 2*basis[0][0] # 2 for hexagonal self
#            height = self.basis[1][1]

#            tmp=layer+(0,slice(2))
#            temp = self.coordinates_to_orbital_position[tmp]
#            for a_x, a_y in zip(*(temp.reshape([self.dimensions[0]*self.dimensions[1],2])+self.com[:2]).T):
#                self.ax.add_patch(Rectangle(
#                    xy=(a_x-width/2, a_y-height/2) ,width=width, height=height,
#                    linewidth=1, color='blue', fill=False))
        
#        # Annotate orbitals:
#        for o, orb in enumerate(self.orbitals):
#            if annotate_orbs==True:
#                self.ax.annotate(self.label_orb[o], self.orbitals[o][:2]+np.array([0.02,0.02]), horizontalalignment='right', verticalalignment='bottom', color=self.color_orb[o])

#        # Draw basis vectors:
#        if show_basis_vectors==True:
#            for i in range(2): # Two basis vectors are plotted in the plane
#                self.ax.annotate(text='', xytext=self.basis[i][:2],xy=self.orbitals[0][:2],  arrowprops=dict(arrowstyle='<-', lw=2))
#                self.ax.annotate(text=self.label_basis[i], xytext=0.5*self.basis[i][:2]+[0.05,0.05],xy=self.orbitals[0][:2])

#        return self.fig, self.ax

#    def _plt_hopping(self, orb_i: int, orb_f: int, cell_hop: tuple, label='', PBC=True, show_only_centre=True):
#        cell_hop = np.array(cell_hop)
#        cells = self.n_cells
#        if (cell_hop==np.zeros([self.n_dim])).all() and PBC==True:
#            if show_only_centre:
#                index=[0]
#            else:
#                index=np.arange(cells)
#            for i in index:
#                index_i = self.coord[i][:self.n_dim]
#                index_f = np.mod(index_i + cell_hop, self.dimensions)
#                index_i = list(index_i)+[orb_i]
#                orb_0 = self.coordinates_to_orbital_position[tuple(index_i)] 
#                index_f = list(index_f)+[orb_f]
#                orb_1 = self.coordinates_to_orbital_position[tuple(index_f)]
#                orb_0=orb_0[:2]
#                orb_1=orb_1[:2]
#                self.ax.annotate(text='', xytext=orb_1, xy=orb_0,  arrowprops={'arrowstyle': '<->', 'ls': 'dashed'})
#                self.ax.annotate(text=label, xytext=0.5*(orb_0+orb_1)+(0,0.01), xy=orb_0)
#        # else:
#        for dimension in range(self.n_dim):
#            temp=np.zeros(self.n_dim)
#            temp[dimension]=1
#            if (temp==cell_hop).all():
#                if show_only_centre:
#                    index=[0]
#                else:
#                    index=np.arange(cells)
#                for i in index: 
#                    index_i = self.coord[i][:self.n_dim]
#                    index_f = np.mod(index_i + cell_hop, self.dimensions)
#                    coord_i = self.coord_cell[i]
#                    edge = bool(coord_i[dimension]>=self.centre[dimension])
#                    if edge and not PBC:
#                        pass
#                    else:
#                        index_i = list(index_i)+[orb_i,ALL]
#                        orb_0 = self.coordinates_to_orbital_position[tuple(index_i)] 
#                        index_f = list(index_f)+[orb_f,ALL]
#                        orb_1 = self.coordinates_to_orbital_position[tuple(index_f)]
#                        orb_0=orb_0[:2]
#                        orb_1=orb_1[:2]
#                        self.ax.annotate(text='', xytext=orb_1, xy=orb_0,  arrowprops={'arrowstyle': '<->', 'ls': 'dashed'})
#                        self.ax.annotate(text=label, xytext=0.5*(orb_0+orb_1)+(0,0.01), xy=orb_0)



    ############################################################
    ################### Greens' function #######################
    ############################################################

class GreensFunction(BogoliubovdeGennes):
    """Greens function
    
    Attributes
    ----------

    Parameters
    ----------
    conf : model = instance of TightBinding or BogoloiubovdeGennes

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

    def __init__(self, model, energy_interval, resolution, k_axes='default'):

        self.lattice_vectors = model.lattice_vectors
        self.name = model.name

        self._model= model._model
        
        if k_axes!='default':
            model.set_k_axes(k_axes)

        self._k_axes = np.copy(model._k_axes)

        self._plot_index=0

        self.kpts = model.kpts

        self.n_spins = model.n_spins
        self._atoms = model._atoms
        self._pieces = model._pieces
        self._atom_dict= model._atom_dict
        self._atom_indices = model._atom_indices
        self._counter= model._counter
        self._glue_edgs = model._glue_edgs
        self.bulk_calculation=model.bulk_calculation

        self.temperature=model.temperature

        self.hoppings=model.hoppings
        self.hopping_labels=model.hopping_labels
        self.impurities=model.impurities
        self.impurity_labels=model.impurity_labels

        self._flattened_kpts = model._flattened_kpts

        self.bulk_calculation = model.bulk_calculation

        self.energy_interval=energy_interval
        self.emax=np.max(energy_interval)
        self.emin=np.min(energy_interval)
        self.resolution=resolution
        self.n_energy=len(self.energy_interval)

        self._Berry_curvature = False
        
        if self._model=='tb':
            dm = np.array([np.multiply(model.eigenvectors,np.conj(model.eigenvectors))])
        else:
            density_matrix = np.einsum('...,...->...',model.eigenvectors,np.conj(model.eigenvectors),optimize=True)
            omegas = np.array(self.energy_interval, dtype=COMPLEX) + 1.0j*self.resolution
            reciprocal = 1/(omegas[:,None]-model.eigenvalues[:]) #assuming T=0
            if self.bulk_calculation:
                self._dos = np.einsum('koe,k...e->k...o', reciprocal, density_matrix,optimize=True)
            else:
                self._dos = np.einsum('oe,...e->...o', reciprocal, density_matrix,optimize=True)

        if self._Berry_curvature:
            self.Berry_curvature()

        def Berry_curvature(self):
            ############################################
            ############# Berry curvature ##############
            ############################################
            psi_x=np.conj(model.eigenvectors)
            psi_y=np.conj(model.eigenvectors)
            [xx,yy]=self.centre
            for x in range(-xx,xx+1):
                for y in range(-yy,yy+1):
                    if x==0:
                        dx=0
                    else:
                        dx=int(x/np.abs(x))
                    if y==0:
                        dy=0
                    else:
                        dy=int(y/np.abs(y))
                    psi_x[:,x,y,:,:,:]=(psi_x[:,x+dx,y,:,:,:]-psi_x[:,x-dx,y,:,:,:])/2
                    psi_y[:,x,y,:,:,:]=(psi_y[:,x,y+dy,:,:,:]-psi_y[:,x,y-dy,:,:,:])/2


            density_matrix_x = np.einsum('...,...->...',model.eigenvectors,psi_x,optimize=True)
            density_matrix_y = np.einsum('...,...->...',model.eigenvectors,psi_y,optimize=True)
            temp = (-1/np.pi)*np.imag([density_matrix_x,density_matrix_y])
            # temp = (-1j/np.pi)*np.array([density_matrix_x,density_matrix_y])
            temp = temp[1]
            temp = np.sum(temp,3)
            temp = np.sum(temp,3)
            temp = temp[...,self.n_dof]
            X=np.arange(-self.centre[0],self.centre[0]+1)*2*np.pi/(self._pieces[0])
            Y=np.arange(-self.centre[1],self.centre[1]+1)*2*np.pi/(self._pieces[1])
            U=temp[0]
            V=temp[1]
            U=np.fft.fftshift(U, axes=1)
            V=np.fft.fftshift(V, axes=1)
            C=np.sqrt(U**2+V**2)
            U=U/C
            V=V/C

            with open(DATA+'temp', 'wb') as f:
                cPickle.dump([X,Y,U,V],f)

            [X,Y,U,V]=np.load(DATA+'temp', allow_pickle=True)

            fig, ax = plt.subplots()
            ax.quiver(X,Y,U.T,V.T)#,color='cyan')
            # ax.set_xticks([-np.pi,0,np.pi])
            # ax.set_yticks([-np.pi,0,np.pi])
            # ax.set_xticklabels([f'$-\pi$',f'$0$',f'$\pi$'])
            # ax.set_yticklabels([f'$-\pi$',f'$0$',f'$\pi$'])

            ax.set_xlim(-0.1,0.1)
            ax.set_ylim(0.25,0.4)
            # C=np.fft.fftshift(C,axes=1)
            # ax.imshow(C.T,cmap='inferno',origin='lower')
            plt.show()
            
            exit()
            ############################################

    def _query_dos(self):
        if type(self._dos)==type(None):
            raise ValueError("Density of states doesn't exist. Please run self.calculate_greens_function first!")
        else: 
            pass

#     def calculate_greens_function(self,energy_interval,resolution,ft_axes=None):

#         return self._dos

    def density_of_states(self, sites='resolved', atom='integrated', orbital='integrated', spin='integrated', energy='resolved', anomalous=False):
        """If site=None: trace sites, else local density of states
If atom=None: trace atoms, else atom resolved
If orbital=None: trace orbitals, else orbital resolved
If spin=None: trace spin, else spin polarised"""
        self._query_dos()

        temp = (-1/np.pi)*np.imag(self._dos)

        if anomalous:
            temp=temp[1]
        else:
            temp=temp[0]

        dim=self.n_dimensions
        if spin=='integrated':
            temp = np.sum(temp,-2)
        elif spin=='resolved':
            pass
        else:
            temp = temp[...,self.spin(spin),:]
            # temp = temp[...,self.spin(spin),:,:]
        if atom=='integrated' and orbital=='integrated':
            temp = np.sum(temp,dim)
        elif atom=='resolved' and orbital=='resolved':
            pass
        elif atom=='integrated':
            indices=[]
            for atom in self._atom_indices:
                indices.append(self.index(atom,orbital))
            temp = temp[..., indices,:]
            # temp = temp[..., indices]
            temp = np.sum(temp,dim)
        elif orbital=='integrated':
            indices=[]
            for orbital in self.atom(atom)._orbital_indices:
                indices.append(self.index(atom,orbital))
            temp = temp[..., indices,:]
            # temp = temp[..., indices]
            temp = np.sum(temp,dim)
        else:
            index=(self.index(atom,orbital))
            temp = temp[..., index,:]
        if sites=='resolved':
            pass
        elif sites=='integrated':
            for i in range(self.n_dimensions):
                temp = np.sum(temp,0)
        else:
            tmp=[]
            for site in sites:
                for coord in site:
                    tmp.append(temp[coord])
            temp=np.sum(tmp,0)
        if energy=='integrated':
            temp = np.sum(temp,-1)
        elif energy=='resolved':
            pass
        else:
            index = FindNearestValueOfArray(self.energy_interval,energy)
            temp = temp[..., index]
        return temp

    def integrated_density_of_states(self, sites='resolved', atom='integrated', orbital='integrated', spin='integrated', anomalous=False):
        return self.density_of_states(sites=sites, atom=atom, orbital=orbital, spin=spin, energy='integrated', anomalous=anomalous)

    def local_density_of_states(self, sites='resolved', atom='integrated', orbital='integrated', energy='resolved', anomalous=False):
        return self.density_of_states(sites=sites, atom=atom, orbital=orbital, spin='integrated', energy=energy, anomalous=anomalous)

    def spin_polarised_local_density_of_states(self, sites='resolved', atom='integrated', spin='resolved', energy='resolved', anomalous=False):
        return self.density_of_states(sites=sites, atom=atom, orbital='integrated', spin=spin, energy=energy, anomalous=anomalous)

    def local_density_of_states(self, sites='resolved', atom='integrated', energy='resolved', anomalous=False):
        return self.density_of_states(sites=sites, atom=atom, orbital='integrated', spin='integrated', energy=energy, anomalous=anomalous)

    def spin_polarised_local_density_of_states(self, sites='resolved', atom='integrated', spin='resolved', energy='resolved', anomalous=False):
        return self.density_of_states(sites=sites, atom=atom, orbital='integrated', spin=spin, energy=energy, anomalous=anomalous)

    def staggered_density(self, atom_i, atom_f, energy=None):
        A=self.local_density_of_states(energy=energy, atom=atom_i, orbital=None, spin=None)
        B=self.local_density_of_states(energy=energy, atom=atom_f, orbital=None, spin=None)
        return A-B

    def mean_abs_staggered_density(self, atom_i, atom_f, energy=None):
        den=self.staggered_density(atom_i,atom_f,energy)
        den=np.abs(den)
        for i in range(self.n_dimensions):
            den=np.sum(den,axis=0)
        den=den/self.n_cells
        return den

    def IPR_abs_staggered_density(self, atom_i, atom_f, energy=None):
        stag_density=self.staggered_density(atom_i=atom_i,atom_f=atom_f,energy=energy)
        ipr=np.sum(np.abs(stag_density))**2/np.sum(np.square(stag_density))
        return ipr

    def magnetism(self, energy=None, atom=None, orbital=None):
        up=self.local_density_of_states(energy, atom, orbital, spin=self.spin('up'))
        dn=self.local_density_of_states(energy, atom, orbital, spin=self.spin('down'))
        return up-dn

    def mean_abs_magnetism(self, energy=None, atom=None, orbital=None):
        magnetism=self.magnetism(energy, atom, orbital)
        magnetism=np.abs(magnetism)
        for i in range(self.n_dimensions):
            magnetism=np.sum(magnetism,axis=0)
        magnetism=magnetism/self.n_cells
        return magnetism

    def spectrum(self, locations, orbital, spins=None):
        self._query_dos()
        temp = self._dos
        if type(spin)==type(None):
            temp = np.sum(temp,-2)
        else:
            temp = temp[...,spin,:]
        temp = temp[..., orbital,:]
        location=locations[0]
        temp = np.array([temp[(*location, ...)] for location in locations])
        return temp

    #############################
    ####### Check these:  #######
    #############################
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

    def plot_spectrum(self, ax, energy='resolved', axes=['resolved','integrated'], atom='integrated', xmin='default', xmax='default', omega_min='default', omega_max='default', vmin='default', vmax='default',label=''):

        from matplotlib.ticker import FuncFormatter, MultipleLocator
        ldos=self.local_density_of_states(energy='resolved', atom=atom)

        axis_index=axes.index('resolved')
        rlabels=[r'$x/|b_0|$',r'$y/|b_1|$',r'$z/|b_2|$']
        klabels=[r'$k_x|b_0|$',r'$k_y/|b_1|$',r'$k_z/|b_2|$']
        rrlabels=[r'$x$',r'$y$',r'$z$']
        kklabels=[r'$k_x$',r'$k_y$',r'$k_z$']
        axes_labels=[]
        coords=[]
        integrated=[]
        k=0
        for i in range(self.n_dimensions):
            if self._k_axes[i]:
                axes_labels.append(klabels[i])
                if type(axes[i])!=str:
                    kklabels[i]=kklabels[i]+f'$={(axes[i])}$'
                coords.append(kklabels[i])
            else:
                axes_labels.append(rlabels[i])
                if type(axes[i])!=str:
                    rrlabels[i]=rrlabels[i]+f'$={(axes[i])}$'
                coords.append(rrlabels[i])
            if axes[i]=='resolved':
                if self._k_axes[i]:
                    xxmin,xxmax=-np.pi,np.pi
                    if xmin=='default' and xmax=='default':
                        xmin,xmax=xxmin,xxmax
                        j=0
                    else:
                        j=1
                else:
                    xxmin,xxmax=-self.centre[1-i],self.centre[1-i]
                    if xmin=='default' and xmax=='default':
                        xmin,xmax=xxmin,xxmax
                    j=1
            elif axes[i]=='integrated':
                ldos=np.sum(ldos,i-k)
                k+=1
                if self._k_axes[i]:
                    integrated.append(kklabels[i])
                else:
                    integrated.append(rrlabels[i])
            else:
                k_index = FindNearestValueOfArray(np.linspace(-np.pi,np.pi,self._pieces[i]),axes[i])-self.centre[i]
                if i-k==0:
                    ldos=ldos[k_index]
                elif i-k==1:
                    ldos=ldos[:,k_index]
                elif i-k==2:
                    ldos=ldos[:,:,k_index]

        if omega_min=='default':
            omega_min=self.emin
        if omega_max=='default':
            omega_max=self.emax

        ldos=np.fft.fftshift(ldos.T, axes=1)

        if energy!='resolved':
            energy_index = FindNearestValueOfArray(self.energy_interval,energy)
            ldos=ldos[energy_index]

        if vmin=='default':
            vmin=np.min(ldos)
        if vmax=='default':
            vmax=np.max(ldos)

        ax.set_xlim(xmin,xmax)
        ax.set_xticks([xmin,0,xmax])
        if j==0:
            ax.set_xticklabels([f'$-\pi$',f'$0$',f'$\pi$'])
        elif j==1:
            ax.set_xticklabels([f'${xmin}$',f'$0$',f'${xmax}$'])
        ax.set_xlabel(axes_labels[axis_index])
        ax.set_ylabel('$\omega$')

        title=''
        coords=', '.join(coords)
        if len(integrated)>0:
            integrated=', '.join(integrated)
            title=r'$\int\text{{d}}$'+integrated
        title=title+r'$-\frac{{1}}{{\pi}}\Im\hat{{\mathcal{{G}}}}(\omega-i\epsilon,$ '+coords+'$)$'
        colors=['r^-','bo-','g*-']
        if energy=='resolved':
            title=title+f'$|_{{\epsilon={self.resolution}}}$'
            ax.set_ylim(omega_min,omega_max)
            ax.set_yticks([omega_min,0,omega_max])
            extent=[xmin,xmax,self.emin,self.emax]
            ax.imshow(ldos, extent=extent, origin='lower', vmin=vmin, vmax=vmax,aspect='auto', interpolation=None)
        else:
            title=title+f'$|_{{\omega={energy}, \epsilon={self.resolution}}}$'
            ax.set_ylim(vmin,vmax)
            ax.set_yticks([vmin,vmax])
            xs=np.linspace(xxmin,xxmax,self._pieces[axis_index])
            ax.plot(xs,ldos,colors[self._plot_index],label=label)
            self._plot_index+=1
            ax.set_ylabel('Spectral density')
        ax.set_title(title)

        return ax

    def plot_ldos(self, ax, energy='resolved', axes=['resolved','resolved'], atom='integrated', xmin='default', xmax='default', ymin='default', ymax='default', vmin='default', vmax='default',label=''):

        ################################
        from matplotlib.ticker import FuncFormatter, MultipleLocator
        ldos=self.local_density_of_states(energy=energy, atom=atom)
        ldos=np.fft.fftshift(ldos.T)
        
        if 'integrated' in axes:
            axis_index=axes.index('integrated')
        rlabels=[r'$x/|b_0|$',r'$y/|b_1|$',r'$z/|b_2|$']
        klabels=[r'$k_x|b_0|$',r'$k_y/|b_1|$',r'$k_z/|b_2|$']
        rrlabels=[r'$x$',r'$y$',r'$z$']
        kklabels=[r'$k_x$',r'$k_y$',r'$k_z$']
        axes_labels=[]
        coords=[]
        integrated=[]
        k=0
        l=0
        for i in range(self.n_dimensions):
            if self._k_axes[i]:
                axes_labels.append(klabels[i])
                if type(axes[i])!=str:
                    kklabels[i]=kklabels[i]+f'$={(axes[i])}$'
                coords.append(kklabels[i])
            else:
                axes_labels.append(rlabels[i])
                if type(axes[i])!=str:
                    rrlabels[i]=rrlabels[i]+f'$={(axes[i])}$'
                coords.append(rrlabels[i])
            if axes[i]=='resolved':
                if l==0:
                    if self._k_axes[i]:
                        xxmin,xxmax=-np.pi,np.pi
                        if xmin=='default' and xmax=='default':
                            xmin,xmax=xxmin,xxmax
                            jx=0
                        else:
                            jx=1
                    else:
                        xxmin,xxmax=-self.centre[1-i],self.centre[1-i]
                        if xmin=='default' and xmax=='default':
                            xmin,xmax=xxmin,xxmax
                        jx=1
                    l+=1
                if l==1:
                    if self._k_axes[i]:
                        yymin,yymax=-np.pi,np.pi
                        if ymin=='default' and ymax=='default':
                            ymin,ymax=yymin,yymax
                            jy=0
                        else:
                            jy=1
                    else:
                        yymin,yymax=-self.centre[1-i],self.centre[1-i]
                        if ymin=='default' and ymax=='default':
                            ymin,ymax=yymin,yymax
                        jy=1
            elif axes[i]=='integrated':
                ldos=np.sum(ldos,i-k)
                k+=1
                if self._k_axes[i]:
                    integrated.append(kklabels[i])
                else:
                    integrated.append(rrlabels[i])
            else:
                k_index = FindNearestValueOfArray(np.linspace(-np.pi,np.pi,self._pieces[i]),axes[i])-self.centre[i]
                if i-k==0:
                    ldos=ldos[k_index]
                elif i-k==1:
                    ldos=ldos[:,k_index]
                elif i-k==2:
                    ldos=ldos[:,:,k_index]

        if vmin=='default':
            vmin=np.min(ldos)
        if vmax=='default':
            vmax=np.max(ldos)

        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_xticks([xmin,0,xmax])
        ax.set_yticks([ymin,0,ymax])
        if jx==0:
            ax.set_xticklabels([f'$-\pi$',f'$0$',f'$\pi$'])
        elif jx==1:
            ax.set_xticklabels([f'${xmin}$',f'$0$',f'${xmax}$'])
        if jy==0:
            ax.set_yticklabels([f'$-\pi$',f'$0$',f'$\pi$'])
        elif jy==1:
            ax.set_yticklabels([f'${ymin}$',f'$0$',f'${ymax}$'])

        ax.set_xlabel(axes_labels[0])
        ax.set_ylabel(axes_labels[1])
        self.xlabel=axes_labels[0]
        self.ylabel=axes_labels[1]

        title=''
        coords=', '.join(coords)
        if len(integrated)>0:
            integrated=', '.join(integrated)
            title=r'$\int\text{{d}}$'+integrated
        title=title+r'$-\frac{{1}}{{\pi}}\Im\hat{{\mathcal{{G}}}}(\omega-i\epsilon,$ '+coords+'$)$'
        colors=['r^-','bo-','g*-']
        extent=[xmin,xmax,ymin,ymax]
        if energy=='resolved':
            title=title+f'$|_{{\epsilon={self.resolution}}}$'
        else:
            title=title+f'$|_{{\omega={energy}, \epsilon={self.resolution}}}$'
        ax.set_title(title)
        self.title=title
        ################################

        ax.imshow(ldos, extent=extent, origin='lower', vmin=vmin, vmax=vmax, aspect='auto', interpolation=None)
        return ax
