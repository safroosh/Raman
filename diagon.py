"""
#######################################################
## Author: Sajjad Afroosheh                          ##
## email: sajjada@bgsu.edu                           ##
## Zayak's research group                            ##
## http://physics.bgsu.edu/%7Eazayak/Welcome.html    ##
## Bowling Green State University                    ##
#######################################################

"""


"""
    Input: file of forces due to positive displacement
           file of forces due to negative displacment
           ---required format: each line of file corresponds to single displacement,
           no punctuation, each force separated by space
           Number of total atoms, number of each type of atom 
           mass of each atom (given in order that appears in scf file)
           
    Function: uses input forces to determine vibrational frequencies
    
    Output: eigenvectors and eigenvalues of dynamic matrix created from input data
        vibrational frequencies       
"""

"""
    Vibrations of the structures were computed by the procedure of Postnikov:
        
            DOI:	10.1103/PhysRevB.71.115206
"""

"""
    In actual macro, the structural data would be read from
    structutral input file. The data in the file is considered
    as intellectual proprerties, hence, a dummy input considered
    in this demo.
"""

import sys
#import os
import periodictable as PT
import numpy as np
#import cmath
#import re
#import sys
import math
from numpy import linalg as LA
#caution: atomic displacement of 0.2 A assumed as well as inital force units of ry/au

def DynMtx():
    np.set_printoptions(threshold=sys.maxsize) #avoid truncation on future printing
    #   pos, neg, struc = readInp()
    neg = np.array(read_struc('forces_minus.dat'))
            #forces_plus = read_forces('forces_plus')
    pos = np.array(read_struc('forces_plus.dat'))
    struc = read_struc('siesta.cor') #use siesta code to read atoms


    """forces will be loaded in having units of ry/au
    convert ry/au to eV/A (Ry/au = 25.711 eV/A)"""
    # Conversion factors and displacement in meters
 #   ry_to_ev = 25.711  # 1 Ry = 25.711 eV
    ev_to_N = 1.60217e-9  
    angstrom_to_m = 1e-10
    displacement = 0.2 * angstrom_to_m

    # Convert forces from Ry/au to N
    pos *= ev_to_N 
    neg *= ev_to_N
    """ Gather structural information
    #requests for total number of atoms in the system
    #variables initialized for loop"""
    natoms = len(struc) #int(struc[:,0].sum())
    dimension = 3 * natoms        #dimension of system (xyz-coordinates for each atom)
#    masses = [struc[i,1] for i in range(len(struc[:,0])) for _ in range(int(struc[i,0]))]
    avg = (pos - neg) / (2 * displacement) #* 1.889725989
    
    D = construct_dynamic_matrix(avg, dimension, struc)
    
    freq, eigV = calculate_frequencies(D)
    print(freq/10000000000)
    
    write_configurations(struc, eigV, freq)
    write_frequencies(freq)
    

    # Output frequencies and modes
 #   writeOut(freq, eigV)


#==================================================================================
def construct_dynamic_matrix(avg, dimension, struc):
    """
    Constructs the dynamic matrix based on averaged forces and atom masses.
    
    begins construction of dynamic matrix
    averages two force matrices, symmetrizes, divides by sqrt of product of masses
    converts to inverse centimeters
    
    """
    D = np.zeros((dimension, dimension), dtype=float)
    
    for x in range(dimension):
        for y in range(x, dimension):
            # Determine which atoms are in consideration
            massa, massb = getattr(PT,struc[x // 3][-1]).mass, getattr(PT,struc[y // 3][-1]).mass #determines which molecules are in consideration
            # Calculate symmetric element of the dynamic matrix
            d = -0.5 * (avg[x, y] + avg[y, x])  # eV/A^2
            
            # Calculate mass product and adjust dynamic matrix element
            mass_product = math.sqrt(massa * massb) * 1.66e-24 #divides out masses
            D[x, y] = D[y, x] = d / mass_product  # eV/m^2/kg #sets to corresponding dynamic matrix element in #eV/gA^2

    return D
#==================================================================================
def calculate_frequencies(D):
    """
    Calculates vibrational frequencies from the dynamic matrix.
    """
    # Calculate eigenvalues and eigenvectors. # Diagonalize dynamic matrix to find eigenvalues and eigenvectors
    eig, eigV = LA.eigh(D)
    
    # Skip the first 6 eigenvalues and eigenvectors (translational and rotational modes)
    eig, eigV = eig[6:], eigV[6:,:]
    
    # Convert eigenvalues to frequencies in THz
    freq = np.sqrt(eig) * 33.35641 / (2 * np.pi)
    
    
    return freq, eigV

        
#==================================================================================

def convert(item):
    try:
        # Try converting to float
        return float(item)
    except ValueError:
        # Return as string if conversion fails
        return item
    
def read_struc(file_name):
    
    with open(file_name, 'r') as file:
        content = []
        for line in file:
        # Apply conversion to each item in the line
            content.append([convert(number) for number in line.split()])
    return content


#==========================================================================================

def write_configurations(Atoms, eigV, freq):
    """
    Writes configuration files for positive and negative displacements
    based on eigenvalues and original coordinates.
    
    Atoms: list of strings containing atom coordinates and types.
    eigV: Matrix of eigenvectors.
    freq: List of frequencies.
    """

    coords = []
    atom_types = []
    atom_symbols = []

# Iterate once over Atoms to extract all necessary information
    for atom in Atoms:
        coords.append(atom[:3])
        atom_types.append(atom[3])
        atom_symbols.append(atom[-1])

# Convert coords to a NumPy array in one step
    coords = np.array(coords, dtype=float)
 #   formatted_eigV = np.array(eigV).reshape(-1, len(coords), 3)[6:]
    
    for i in range(len(coords)):

        formatted_eigV = eigV[i].reshape((len(coords), 3))
        new_pos_m = coords - formatted_eigV #formatted_eigV[i]
        new_pos_p = coords + formatted_eigV #formatted_eigV[i]

            
        with open(f'm_{freq[i]:.5f}.xyz', 'w') as f,  open(f'Conf_m{i+7}.dat', 'w') as fm, open(f'Conf_p{i+7}.dat', 'w') as fp:
    
            f.write(f'{len(coords)}\n\n')
    
            lines = []
            Confm =[]
            Confp =[]
            
            for coord, eigV1, content_m, content_p, symbol, atype in zip(coords, formatted_eigV, new_pos_m, new_pos_p, atom_symbols, atom_types):
                line = f"{coord[0]:7.5f} {coord[1]:7.5f} {coord[2]:7.5f}  {eigV1[0]:<13.5e} {eigV1[1]:<13.5e} {eigV1[2]:<13.5e}   {symbol}"
                confm = f"{content_m[0]:7.5f} {content_m[1]:7.5f} {content_m[2]:>7.5f} {int(atype)} {symbol}"
                confp = f"{content_p[0]:7.5f} {content_p[1]:7.5f} {content_p[2]:>7.5f} {int(atype)} {symbol}"
                
                lines.append(line)
                Confm.append(confm)
                Confp.append(confp)

            f.write('\n'.join(lines))
            fm.write('\n'.join(Confm))
            fp.write('\n'.join(Confp))

   
    write_mod_nmd(coords,  atom_symbols, eigV)
            
                

def write_frequencies(freq):
    """
    Writes frequencies to a file.
    """
    with open('freq.dat', 'w') as f:
        for value in freq:
            f.write(f"{value:.5f}\n")


def write_mod_nmd(coords,  atom_symbols, eigV):
    """
    Writes a single file containing all modes for all atoms.
    """
    with open('mod.nmd', 'w') as f:
        f.write('atomnames ' + ' '.join(atom_symbols) + '\n')
        f.write('coordinates ')
        for i in range(len(coords)):
            f.write(' '.join(map(str, coords[i])))
        for i in range(len(eigV)):
            f.write(f'\nmode {i+7} ')
            f.write(' '.join(map(str, eigV[i])))


if __name__ == "__main__":

    DynMtx()
    





