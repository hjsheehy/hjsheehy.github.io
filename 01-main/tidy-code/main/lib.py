import sys
import os
from importlib.util import spec_from_loader, module_from_spec
from importlib.machinery import SourceFileLoader
import numpy as np
import scipy
from scipy import linalg as la
import scipy.sparse as sp
from math import floor, ceil
import glob
import time
import zipfile as zp
import gc
'''A library for multiorbital non-unitary spin-triplet superconductivity on a lattice
simulation.'''
###################################################################
########################## Scripting ##############################
###################################################################
def import_path(path):
    '''e.g. path = 'subfolder/file.conf' 
    Imports entire script with arbitrary file extension'''
    module_name = os.path.basename(path)#.replace('-', '_')
    spec = spec_from_loader(
        module_name,
        SourceFileLoader(module_name, path)
    )
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[module_name] = module
    return module
#
def import_all(module):
    '''Equivalent to 'from module import *' for modules with extensions
    other than .py'''
    module_dict = module.__dict__
    try:
        to_import = module.__all__
    except AttributeError:
        to_import = [name for name in module_dict if not name.startswith('_')]
    return module_dict, to_import
###################################################################
########################### General ###############################
###################################################################
Pauli_z = np.array([[1,0],[0,-1]])
Pauli_y = np.array([[0,-1.0j],[1.0j,0]])
Pauli_x = np.array([[0,1],[1,0]])
Pauli_plus = (Pauli_x+1.0j*Pauli_y)/2
Pauli_minus = (Pauli_x-1.0j*Pauli_y)/2
Pauli_vec = np.dstack([Pauli_x, Pauli_y, Pauli_z])
#
def Dagger(array):
    array = np.conj(array.T)
    return array
#
def Trace(gr):
    '''Takes the (diagonalised) Green's function GR and returns the Trace over
    the electron part and spin and orbital components'''
    # hole=np.array([[0,0],[0,1]]) #hole part
    electron=np.array([[1,0],[0,0]]) #electron part
    spins=orbitals=np.eye(2)
    index=np.diag(np.kron(electron, np.kron(spins, orbitals)))
    return np.inner(index, gr)
#
def T_rev(v, dof):
    temp = np.zeros([2*dof,2*dof],dtype=complex)
    temp[0:dof,:] = -np.conj(v[dof:2*dof,:])
    temp[dof:2*dof,:] = np.conj(v[0:dof,:])
    return temp
###################################################################
######################### Hamiltonian #############################
###################################################################
def NN(site: tuple, n: tuple, nn_vectors: "list of tuples") -> "list of tuples":
    '''Forward nearest-neighbours of an element of the lattice group'''
    nn=[]
    dimensions=len(n)
    for i in range(dimensions):
        nn.append(site)
    return nn
#
def H0(lattice_basis: "list of tuples",
        site_tensors: "list of kron(sites,spin,orbital)",
        hopping: "list of lists of tuples",
        hopping_tensors: "list of lists of kron(sites,spin,orbital)",
        PBC: "list of lists of Boolean",
        impurity_locations: "list of tuples",
        impurity_tensors: "list of kron(spin,orbital)",
        dimensions: tuple,
        n_spins: int,
        n_orbitals: int) -> "complex128 2D array":
    '''Multiorbital spin-triplet tight-binding Hamiltonian with impurities
    \nThe hopping_vectors are of the form [[intracell],[inter-nearest-cell],[inter-next-nearest cell],...]'''
    n_dim = len(lattice_basis)
    n_sites = len(site_tensors)
    n_cells = np.prod(dimensions)
    n_atoms = n_sites * n_cells
    dof = n_atoms * n_spins * n_orbitals
    onsite = np.kron(np.eye(n_atoms),onsite_tensor)
    n_nn = sum(len(level) for level in hopping_vectors)
    # hopping
    graph = np.zeros([n_nn, n_sites, n_sites])
    hop = np.zeros([n_nn, dof, dof]) # spin-orbit
    intracell_hopping = hopping[0]
    n_intracell=len(intracell_hopping)
    intercell_hopping = hopping[1]
    n_intercell=len(intercell_hopping)
    for i in range(n_intracell):
        site_i=intracell_hopping[i][0]
        site_f=intracell_hopping[i][1]
        for j in range(n_cells):
            index_i=n_cells*site_i+j
            index_f=n_cells*site_j+j
            hop[i,index_i,index_f]=1

    for nn in range(n_nn):
        for cell in range(n_cells):
            site_A=hopping[
            r = Z_Zn(site, n_dim)
            nn = r + nn_vector        
            nn_int = Zn_Z(nn, n_dim)
            graph[dim,site,nn_int]=1
            if (r[dim]==(n_dim[dim]-1)) and not(PBC[i][j]):
                r0=np.copy(r)
                r0[dim]=0
                site0=Zn_Z(r0,n_dim)
                graph[dim,site,site0]=0
            # include spin-orbit to the hopping:
            hop[dim] = np.kron(graph[dim,:,:],hopping_tensors[i][j])
            j+=1
        i+=1
    del(graph)
    hop = np.sum(hop, axis=0) 
    hop += Dagger(hop) #to account for possibility of complex hopping
    # Impurities
    impurities=np.zeros([n_sites,n_sites])
    for impurity in impurity_locations:
        loc=Zn_Z(centre(impurity,n), n)
        impurities[loc,loc] = 1
    impurities=np.kron(impurities,impurity_tensor)
    return onsite + hop + impurities
#
def Mean_field_all(SC_tensor, n_sites):
    '''Spin, orbital pairing mean-field'''
    return np.kron(np.eye(n_sites), SC_tensor)
#
def Mean_field_layer(SC_tensor,layer, n_sites):
    '''Spin, orbital pairing mean-field
    layer :: int : refers to layer in z-direction'''
    mean_field = np.zeros([n_sites])
    for i in range(n_sites):
        site = Z_Zns(i, n_x, n_y, n_z)
        if site[2]==layer:
            mean_field[i]=1
    return np.kron(np.diag(mean_field), SC_tensor)
#
def Mean_field_sites(SC_tensor, sites, n_sites):
    '''Spin, orbital pairing mean-field
    sites :: list'''
    mean_field = np.zeros([n_sites])
    for i in range(n_sites):
        mean_field[i]
    return np.kron(np.eye(n_sites), SC_tensor)
#
def H(H0, mean_field):
    '''BdG Hamiltonian to diagonalise the BCS Hamiltonian in the mean-field 
    limit'''
    # return (sp.kron(Pauli_z, H0)
            # + sp.kron(Pauli_plus, Mean_field)
            # + sp.kron(Pauli_minus, np.conj(Mean_field) )
            # ).todense()
    temp = np.kron(Pauli_plus, mean_field)
    temp += np.kron(Pauli_minus, np.conj(mean_field))
    temp += np.kron(Pauli_z, H0)
    return temp
# New
def H_BdG(h0, Hartree, Fock, Gorkov, index,
            n_x,n_y,n_sites,n_spins,n_orbitals,dof):
    hartree = np.zeros([dof,dof])
    fock = np.zeros([dof,dof],dtype='complex')
    gorkov = np.zeros([dof,dof],dtype='complex')
    part=0
    for x in range(n_x):
        for y in range(n_y):
            for s in range(n_spins):
                for m in range(n_orbitals):
                    index0=index[part, x, y, s, m]
                    for t in range(n_spins):
                        for n in range(n_orbitals):
                            index1=index[part, x, y, t, n]
                            fock[index0,index1] = Fock[x,y,s,t,m,n]
                            gorkov[index0,index1] = Gorkov[x,y,s,t,m,n]
                    index1=index0
                    hartree[index0,index1] = Hartree[x,y,s,m]
    temp = np.kron(Pauli_minus, -np.conj(gorkov) )
    temp +=np.kron(Pauli_plus, gorkov)
    temp +=np.kron(np.array([[0,0],[0,1]]), -np.conj(h0+hartree+fock))
    temp +=np.kron(np.array([[1,0],[0,0]]), h0+hartree+fock)
    return temp
#
def Zn_Z(site : tuple, dimensions : tuple) -> int:
    '''Mapping Zn->Z from the lattice group to the integers'''
    dimension=len(dimensions)
    temp=np.zeros([dimension])
    for i in range(dimension):
        temp[i]=site[i]%dimensions[i]
        for j in range(i):
            temp[i]*=dimensions[j]
    index = int(np.sum(temp))
    return index 
#
def Z_Zn(index : int, dimensions : tuple) -> tuple:
    '''Mapping Z->Zn from the integers to the lattice group'''
    dimension=len(dimensions)
    site=np.zeros([dimension],int)
    for i in range(dimension):
        temp=1
        for j in range(i):
            temp*=dimensions[j]
        site[i]=int(floor(index / temp) %dimensions[i])
    return site
#
def Index(dimensions, n_sites, n_spins, n_orbitals):
    '''Returns an array of the indexing'''
    n_cells = np.prod(dimensions)
    n_atoms = n_cells*n_sites
    dof = n_atoms*n_spins*n_orbitals
    return np.array(
            [[[[[[[Zn_Z([part,x,y,z,a,s,m],
                [2,*dimensions,n_sites,n_spins,n_orbitals])
            for m in range(n_orbitals)]
            for s in range(n_spins)]
            for a in range(n_sites)]
            for z in range(dimensions[2])]
            for y in range(dimensions[1])]
            for x in range(dimensions[0])]
            for part in range(2)], dtype=int)

#
def centre(r,n):
    '''Converts r in central coordinate system to a pair of coordinates
    in the Hamiltonian matrix e.g. r=[0,0] for an impurity at the centre.'''
    return [(r[0]+int(n[0]/2))%n[0], (r[1]+int(n[1]/2))%n[1],r[2]]
###################################################################
############# Density matrix and Greens functions #################
###################################################################
def Density_normal(v, dof):
    '''Electron-electron part. Since real, we change dtype to float.'''
    return (np.real(np.multiply(v,np.conj(v)))).astype('float')
#
# def Density_holes(v, dof):
    # '''Hole-hole part. Since real, we change dtype to float.'''
    # Tv=T_rev(v, dof)
    # return (np.real(np.multiply(Tv,np.conj(Tv)))).astype('float')
#
def Density_anomalous(v, dof):
    '''Magnitude of the SC order parameter. Since real, we change dtype to float.'''
    return (np.abs(np.multiply(v,np.conj(T_rev(v, dof))))).astype('float')
#
# def Density_matrix(density,n_x,n_y,index,n_sites,n_spins,n_orbitals):
    # '''5-dimensional data set: [x, y, spin, orbital, n]'''
    # part = 0 #electron-electron
    # dm = np.array(
            # [[[[density[index[part,x, y, spin, orbital]]
            # for orbital in range(n_orbitals)] for spin in range(n_spins)]
            # for y in range(n_y)] for x in range(n_x)])
    # return dm
############################## NEW: ##############################
def DM(v, n, n_spins, n_orbitals, index):
    n_sites = np.prod(n)
    dof = n_sites*n_spins*n_orbitals

    part = 0 #electron-electron
    # dof=n_sites*n_spins*n_orbitals
    # dm = np.array(
            # [[[[[[(v[indexing(Zn_Z(x, y, n_x, n_y), tau,   n, part, dof, n_sites,n_spins,n_orbitals),:]*
           # np.conj(v[indexing(Zn_Z(x, y, n_x, n_y), sigma, m, part, dof, n_sites,n_spins,n_orbitals),:]))
            # for n in range(n_orbitals)] for m in range(n_orbitals)]
            # for tau in range(n_spins)] for sigma in range(n_spins)]
            # for y in range(n_y)] for x in range(n_x)])
    temp = np.array(
            [[[[[v[index[part,x, y, z, sigma, m],:]
            for m in range(n_orbitals)]
            for sigma in range(n_spins)]
            for z in range(n[2])] for y in range(n[1])] for x in range(n[0])])
    dm = np.einsum('xyzsme,xyztne->xyzstmne',temp,np.conj(temp))
    return dm
#
def Normal_Fermi(dm, dof, w, T):
    if T==0:
        return np.sum(dm[:,:,:,:,:,:,:dof],-1)
    else:
        f = 1/(np.exp(w/T)+1)
        return np.einsum('xyzstmne,e->xyzstmn', dm, f)
#
def ADM(v, n, index):
    n_x, n_y, n_z = N_x(n), N_y(n), N_z(n)
    n_sites = N_sites(n)
    n_spins = N_spins(n)
    n_orbitals = N_orbitals(n)
    n_sites = N_sites(n)
    dof = Dof(n)

    part = 0 #electron-electron
    # dof=n_sites*n_spins*n_orbitals
    # Tv=T_rev(v, dof)
    # dm = np.array(
            # [[[[[[(v[indexing(Zn_Z(x, y, n_x, n_y), tau, n, part, dof, n_sites,n_spins,n_orbitals),:]*
           # np.conj(Tv)[indexing(Zn_Z(x, y, n_x, n_y), sigma, m, part, dof, n_sites,n_spins,n_orbitals),:])
            # for n in range(n_orbitals)] for m in range(n_orbitals)]
            # for tau in range(n_spins)] for sigma in range(n_spins)]
            # for y in range(n_y)] for x in range(n_x)])
    temp = np.array(
            [[[[[[v[index[part,x, y, z, sigma, m],:]
            for m in range(n_orbitals)]
            for sigma in range(n_spins)]
            for z in range(n_z)] for y in range(n_y)] for x in range(n_x)] for part in range(2)])
    dm = np.einsum('xyzsme,xyztne->xyzstmne',temp[0],temp[1])
    return dm
#
def Anomalous_Fermi(adm, dof, w, T):
    if T==0:
        return -np.sum(adm[:,:,:,:,:,:,:dof],-1)
    else:
        f = 1/(np.exp(w/T)+1)
        return np.einsum('xystmne,e->xystmn', (1/2)*adm, 1-2*f)
# 
def indexed(tensor, n_x,n_y,index,n_sites,n_spins,n_orbitals):
    part=0
    indexed_tensor = np.array(
            [[[[[[tensor[index[part,x, y, tau, n],
                         index[part,x, y, sigma, m]]
            for n in range(n_orbitals)] for m in range(n_orbitals)]
            for tau in range(n_spins)] for sigma in range(n_spins)]
            for y in range(n_y)] for x in range(n_x)])
    return indexed_tensor
#
def HartreeFockGorkov(Hubbard_U, dm, adm):
    Hartree=-np.einsum('xystmn,xyttnn->xysm',Hubbard_U,dm)
    Fock=+np.einsum('xystmn,xystmn->xystmn',Hubbard_U,dm)
    Gorkov=-np.einsum('xystmn,xystmn->xystmn',Hubbard_U,adm)
    return np.real_if_close(Hartree), np.real_if_close(Fock), np.real_if_close(Gorkov)
#
# def IndexedHartreeFockGorkov(Hartree_, Fock_, Gorkov_, n_x,n_y,dof,n_sites,n_spins,n_orbitals):
    # Hartree = indexed(Hartree_, n_x,n_y,dof,n_sites,n_spins,n_orbitals)
    # Hartree = np.sum(Hartree, (3,5))
    # Fock = indexed(Fock_, n_x,n_y,dof,n_sites,n_spins,n_orbitals)
    # return Hartree, Fock, Gorkov
#
def Spin_density(Hartree, n_sites, n_orbitals):
    NDown = np.sum(Hartree[:,:,1,:])/(n_sites*n_orbitals)
    NUp = np.sum(Hartree[:,:,0,:])/(n_sites*n_orbitals)
    return NDown, NUp
#
def BCS_Delta(Gorkov, n_sites, n_orbitals):
    Xi = np.sum(Gorkov[:,:,0,1,:,:])/(n_sites*n_orbitals**2) 
    return Xi
##################################################################
def Green_function(omega, eigenvalues, density_matrix):
    '''7-dimensional data set: [x, y, z, spin, spin, orbital, orbital]'''
    return np.einsum('e,xyzstmne->xyzstmn', 1/(omega-eigenvalues), density_matrix)
#
# def Green_function(omega, eigenvalues, density_matrix):
    # '''4-dimensional data set: [x, y, spin, orbital]'''
    # return np.einsum('n,xyosn->xyos', 1/(omega-eigenvalues), density_matrix)
#
def DOS_fixed_energy(omega, eigenvalues, density_matrix):
    '''7-dimensional data set: [x, y, z, spin, spin, orbital, orbital]'''
    green = Green_function(omega, eigenvalues, density_matrix)
    return -(1/np.pi)*np.imag(green)
#   
def DOS(omegas, eigenvalues, density_matrix):
    '''8-dimensional data set: [x, y, z, spin, spin, orbital, orbital, omega]'''
    x = np.array([DOS_fixed_energy(omega, eigenvalues, density_matrix) for omega in omegas])
    dos = np.moveaxis(x,0,-1)
    return dos
###################################################################
######################### Free energy #############################
###################################################################
def BdGGroundState(w):
    '''Finds the ground state of a BdG Hamiltonian, i.e. min positive
    energy. At T=0 this is the free energy.
    w :: eigenvalues'''
    temp_index = np.where(w > 0, w, np.inf).argmin()
    return w[temp_index]
#
# def FreeEnergy():
    
###################################################################
###################### Fourier transform ##########################
###################################################################
def FT(dos):
    '''2D Fourier transform over axes 0 and 1.'''
    f = np.fft.fft2(dos, axes=(0,1), norm='ortho')
    f = np.abs(f)
    fshift = np.fft.fftshift(f, axes=(0,1))
    return fshift
#
def PhaseResolved(sign, dos, omegas, omega):
    '''Input DOS(x,y,spin,orbit,energy)'''
    if len(np.shape(dos))!=5:
        return "Input format :dos: must be DOS(x,y,spin,orbit,energy)"
    temp_index = FindNearestValueOfArray(np.real(omegas), omega)
    f1 = np.fft.fft2(dos[:,:,:,:,temp_index], axes=(0,1), norm='ortho')
    temp_index = FindNearestValueOfArray(np.real(omegas), -omega)
    f2 = np.fft.fft2(dos[:,:,:,:,temp_index], axes=(0,1), norm='ortho')
    f = np.real(f1+sign*f2)
    fshift = np.fft.fftshift(f, axes=(0,1))
    return fshift
#
def PeakPhaseResolved(sign, dos, omegas, omega, dOmega):
    '''Input DOS(x,y,spin,orbit,energy)'''
    if len(np.shape(dos))!=5:
        return "Input format :dos: must be DOS(x,y,spin,orbit,energy)"
    temp_index = FindNearestValueOfArray(np.real(omegas), omega+dOmega)
    f1 = np.fft.fft2(dos[:,:,:,:,temp_index], axes=(0,1), norm='ortho')
    temp_index = FindNearestValueOfArray(np.real(omegas), omega-dOmega)
    f2 = np.fft.fft2(dos[:,:,:,:,temp_index], axes=(0,1), norm='ortho')
    f = np.real(f1+sign*f2)
    fshift = np.fft.fftshift(f, axes=(0,1))
    return fshift

###################################################################
########################### Linecut ###############################
###################################################################
def Straight_line(y, x0, x1):
    x_vals = list(range(x0,x1))
    return [[x, y] for x in x_vals]
#    
def Diagonal_line(centre, edge):
    vals = list(range(centre, edge))
    return [[x, x] for x in vals]
#
def Linecut(DOS_or_IDOS, line):
    '''3-dimensional data set:  DOS(tau, spin, orbital)
    or 4-dimensional data set: IDOS(tau, spin, orbital, omega)
    where tau parameterises the line'''
    data=[]
    for r in line:
        data.append(DOS_or_IDOS[r[0],r[1]])
    return data
###################################################################
##################### Wavefunction basis ##########################
###################################################################
# def Lattice(n_x, n_y, n_sites, basis_lattice, basis_unit_cell):
    # '''input:
    # i :: 2D-array [n_x, n_y] Lattice indices
    # j :: 1D-array [n_sites] Site indices in unit cells
    # basis_lattice :: 2D-array [r_a, r_b]
    # basis_unit_cell :: 2D-array [d_0,..., d_n-1], where n=n_sites'''
    # coordIndex = np.stack([cell for cell in np.ndindex(n_x, n_y)])
    # coordLattice = np.einsum('ij,jk->ik',coordIndex, basis_lattice).reshape(n_x,n_y,2)
    # lattice = np.array([coordLattice + site for site in basis_unit_cell])
    # lattice = np.moveaxis(Lattice, 0, -2)
    # return lattice
#
# def RealLattice(lattice, nx_pts, ny_pts, sigma)
    # '''Gaussian convolution of lattice with Gaussian wavefunctions at
    # each site in each unit cell.    
    # sigma :: scalar
    # Returns: Array of Gaussians [n_x, n_y, n_sites]'''
    
    # norm = 1/np.sqrt(np.pi*sigma**2)
    # return norm*np.exp(-(np.square(np.norm(r-r0)))/sigma**2)

# def GaussianConvolution(indexedArray, real_lattice):
    # '''Convolves an indexed array (e.g. LDOS) with lattice.'''
    # return np.einseum('xyjsme,xyjc->xyjsme', indexedArray, lattice)

###################################################################
########################### Gap map ###############################
###################################################################
def FindNearestValueOfArray(array, value):
    '''Returns index of nearest value in the array'''
    return (np.abs(array - value)).argmin()
#
def FindIndicesOfArray(array, bound1, bound2):
    '''Returns indices of array within boundary'''
    upperBound = max(bound1, bound2)
    lowerBound = min(bound1, bound2)
    return np.where(np.logical_and(array>=lowerBound, array<=upperBound))[0]
#
def CalcSimpleDeriv(y, dx):
    '''y is a 1-d array. returns finite difference vector of length len(v)-1'''
    l=len(y)-1
    y0 = np.delete(y,0,axis=0)
    return np.array([(y0[i]-y[i])/dx for i in range(l)])
#
def FindPeakInSpectrum(inSpectrum, omegas, en1, en2,
    skipresfil, skiprespct, outrangefil, outrangepct):
    nPts = len(inSpectrum)
    dx=abs(omegas[0]-omegas[1])
    
    indices = FindIndicesOfArray(omegas, en1, en2)
    in1,in2 = indices[0], indices[-1]
    
    maxVal = np.amax(inSpectrum[indices])
    maxValIndex = int(np.where(inSpectrum==maxVal)[0])

    peakIndex = -1
    
    posSpectrum = inSpectrum[in1:]
    derivSpectrum = CalcSimpleDeriv(inSpectrum, dx)[in1:]
    decrease = derivSpectrum<0
    
    fiveptsum = decrease \
        + np.append(decrease[1:],[0]) \
        + np.append(decrease[2:],[0,0]) \
        + np.append(decrease[3:],[0,0,0]) \
        + np.append(decrease[4:],[0,0,0,0])
        
    # in0 = FindNearestValueOfArray(fiveptsum, 0)
    
    # indices of points which satisfy 3/5 criterion
    threeoffive = fiveptsum[0:in2-in1]>=3
    count = np.sum(threeoffive)
    
    # if count<=0 then no peak was found between en1 and en2,
    # so return maxValIndex if we suspect the peak is just higher
    # than the given range
    # or return error code -1 if we don’t think the peak is out of range
    if count<=0:
        if outrangefil:
            peakIndex=maxValIndex+in1
        else:
            return peakIndex
    
    threeoffive = np.where(threeoffive)[0]
    # if we are not skipping resonances, then we just return the local max
    # by the first point satisfying the 3of5 criterion
    if not(skipresfil):
        peakIndex = threeoffive[0]
        while ((peakIndex+1)<=(in2-in1)):
            if (posSpectrum[peakIndex]>=posSpectrum[peakIndex+1]):
                break
            else:
                peakIndex = peakIndex + 1
        peakIndex = peakIndex + in1
        
        # if we are skipping resonances, then find the first peak which
        # satisfies the skiprespct criterion
        # if such points exist, then we return the local max around the first one
        # but if such points don’t exist, we consider outrangefil
        peakIndex = threeoffive[0]
        for i in range(count-1):
            tempIndex = threeoffive[i]
            while ((tempIndex+1)<=(in2-in1)):
                if (posSpectrum[tempIndex]>=posSpectrum[tempIndex+1]):
                    break
                else: 
                    tempIndex += 1
            if (posSpectrum[tempIndex]>=posSpectrum[peakIndex]):
                peakIndex = tempIndex
            if (posSpectrum[peakIndex]>=skiprespct/100.0*maxVal):
                break
        if (posSpectrum[peakIndex]>=skiprespct/100.0*maxVal):
            peakIndex = peakIndex + in1
        else:
            # if outrangefil is true then we see if the peakIndex we found
            # (which is the local max by the global max of all pts satisfying
            # the 3of5 criterion) satisfies the outrangepct criterion
            # (which is usually less stringent than the outrangecriterion)
            # if so, then return peakIndex
            # if not, then return maxValIndex
            # if outrangefil is false, then return -1
            if outrangefil:
                if (posSpectrum[peakIndex]>=outrangepct/100.0*maxVal):
                    peakIndex = peakIndex + in1
                else:
                    peakIndex = maxValIndex + in1
            else:
                peakIndex = peakIndex + in1
    return peakIndex
#
def FindGapInSpectrum(inSpectrum, omegas, en1, en2,
        skipresfil, skiprespct, outrangefil, outrangepct):
    return gap
#
def SinglePeakGapMap(spectrum, omegas, en1, en2,
        skipresfil, skiprespct, outrangefil, outrangepct, contiguity, T1, badPixel,
        n_x, n_y, n_spins, n_orbitals):
    gap=np.empty([n_x,n_y,n_spins,n_orbitals])
    # find absolute value of location of the peaks in the spectrum
    # if bad pixel, then set location to -1
    for orbit in range(n_orbitals):
        for spin in range(n_spins):
            for y in range(n_y):
                for x in range(n_x):
                    inSpectrum = spectrum[x,y,spin,spin,orbit,orbit,:]
                    peakIndex = FindPeakInSpectrum(inSpectrum, omegas, en1, en2, 
                                skipresfil, skiprespct, outrangefil, outrangepct)                        
                    if peakIndex==-1:
                        gap[x,y,spin,orbit] = -1
                    else:
                        gap[x,y,spin,orbit] = np.abs(np.real(omegas[peakIndex]))
    # average over bad pixels
    if (contiguity or badPixel):
        for orbit in range(n_orbitals):
            for spin in range(n_spins):
                for y in range(n_y):
                    for x in range(n_x):
                        neighbours = []
                        for i in [-1,1]:
                            temp = gap[(x+i)%n_x,y,spin,orbit]
                            if temp!=-1:
                                neighbours.append(temp)
                            temp = gap[x,(y+i)%n_y,spin,orbit]
                            if temp!=-1:
                                neighbours.append(temp)
                        if len(neighbours)<=0:
                            break
                        av=sum(neighbours)/len(neighbours)
                        a=np.array(neighbours)
                        val=gap[x,y,spin,orbit]
                        if gap[x,y,spin,orbit]==-1:
                            if badPixel == True:
                                gap[x,y,spin,orbit] = av
                        elif contiguity == True:
                            if np.sum(np.where(np.logical_and(a>=val-T1, a<=val+T1))[0])<=2:
                                gap[x,y,spin,orbit] = av
    return gap
#
def GapMap(spectrum, omegas, en1, en2,
        skipresfil, skiprespct, outrangefil, outrangepct, contiguity, T1, badPixel,
        n_x, n_y, n_spins, n_orbitals):
    gapPlus = SinglePeakGapMap(spectrum, omegas, en1, en2,
        skipresfil, skiprespct, outrangefil, outrangepct, contiguity, T1, badPixel,
        n_x, n_y, n_spins, n_orbitals)
    gapMinus = SinglePeakGapMap(spectrum, omegas, -en1, -en2,
        skipresfil, skiprespct, outrangefil, outrangepct, contiguity, T1, badPixel,
        n_x, n_y, n_spins, n_orbitals)
    gapMap=(gapPlus+gapMinus)/2
    for orbit in range(n_orbitals):
        for spin in range(n_spins):
            for y in range(n_y):
                for x in range(n_x):
                    if (gapPlus[x,y,spin,orbit]==-1 or gapMinus[x,y,spin,orbit]==-1):
                        gapMap[x,y,spin,orbit]=-1
    return gapMap
