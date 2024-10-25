import numpy as np



def loadmolecule(file):   # Read xyz file
    f = open(file,'r')
    mol = f.readlines()[2:]  # Skip the two first rows
    atoms = []
    coord = []
    for i in mol:
        atoms.append(i.split()[0])  # Make array of atomic symbols
    coord = np.array([i.split()[1:4] for i in mol], dtype=float)   # Get the coordinates
    return atoms,coord



def mass(atoms):   # Transform atomic symbols array into atomic masses
    atomic_masses = {'H':1.007, 'He':4.002, 'Li':6.938, 'Be':9.012,
                 'B':10.806, 'C':12.009, 'N':14.006, 'O':15.999,
                 'F':18.998, 'Ne':20.180, 'Na':22.989, 'Mg':24.304,    # Dictionary with atomic masses
                 'Al':26.981, 'Si':28.084, 'P':30.973, 'S':32.059,
                 'Cl': 35.446} 
    mol_mass = []
    for j in range(len(atoms)):
        mol_mass.append(atomic_masses[atoms[j]])  # Build the array
    return mol_mass



# Calculate the center of mass
def center_of_mass(masses, coordinates):   
    R=0
    M = sum(m for m in masses)               #Total mass of the system
    for i in range(len(coordinates)):
        R += (masses[i]*coordinates[i]) / M  # Mass-weighted coordinates divided by the total mass
    return R



# Create and diagonalize inertia tensor

def I(newxyz,masses):    # Arguments are coordinates of the centre of mass and the atomic masses array
    mat = np.zeros((3,3),dtype=float)  # Empty matrix
    for i in range(3): # Rows
        for j in range(3): # Columns
            if i == j:    # Target diagonal elements
                for k in range(len(newxyz[:,0])):
                    mat[i,i] += masses[k]*(np.dot(newxyz[k,:],newxyz[k,:])-np.dot(newxyz[k,i],newxyz[k,j]))
                    # I = m*(|r|^2 - alpha^2) with the last being one of the x y or z component
            elif i != j:  # Off-diagonal elements
                for k in range(len(newxyz[:,0])):
                    mat[i,j] += -masses[k]*np.dot(newxyz[k,i],newxyz[k,j]) # I = m*alpha*beta
    return mat



def diag(mat):  # Diagonalize the inertia tensor
    eigval,eigvect = np.linalg.eig(mat)   # Use the numpy function to get eigenvalues and eigenvectors 
    Ia = np.sort(eigval)[2]   # Sort the array with eigenvalues and define the highest as Ia,
    Ib = np.sort(eigval)[1]   # the second largest as Ib and the lowest as Ic
    Ic = np.sort(eigval)[0]    
    D = np.array([[Ic,0,0],[0,Ib,0],[0,0,Ia]]) # Build the diagonal matrix with the eigenvalues in the diagonal
    R = eigvect # The rotation matrix is the eigenvector matrix
    return D,R,Ia,Ib,Ic  # Resturn all the variables



# Recognize symmetric or spherical top

def top(I):
    eigval,eigvect = np.linalg.eig(I) 
    I1 = f'{sorted(eigval)[2]:.2f}'  # Sort and assign eigenvalues
    I2 = f'{sorted(eigval)[1]:.2f}'
    I3 = f'{sorted(eigval)[0]:.2f}'
    # Next are all possible cases regarding the princ. moments of inertia
    if I1 == I2 == I3:       
        print('The molecule is a spherical top')
    elif I1 == I2 != I3:     
        print('The molecule is a prolate symmetrical top')
    elif I1 != I2 == I3:   
        print('The molecule is an oblate symmetrical top')
    elif I1 != I2 != I3:
        print('The molecule is an asymetric top')
    return f''



# Main function that groups the rest

def inertia(molecule): # The argument is the xyz file
    coordinates = loadmolecule(molecule)[1]  # Get coordinates
    atoms = loadmolecule(molecule)[0]        # Get atoms
    atom_mass = mass(atoms)                  # Get masses
    mass_center = center_of_mass(atom_mass,coordinates)   # Compute center of mass
    new_coordinates = coordinates - mass_center           # Translate molecule to the centre of mass
    tensor = I(new_coordinates,atom_mass)       # Build inertia tensor
    diagonal = diag(tensor)[0]                  # Diagonalize tensor and return diag matrix
    I_A = (diag(tensor)[2]*10**(-16))/(6.022*10**23)
    I_B = (diag(tensor)[3]*10**(-16))/(6.022*10**23)  # Get eigenvalues in g*cm^2 units
    I_C = (diag(tensor)[4]*10**(-16))/(6.022*10**23)
    top(tensor)                                  # Classify molecule
    print(f'I_A = {I_A} gcm^2,\nI_B = {I_B} gcm^2, \nI_C = {I_C}gcm^2')  
    rotation = diag(tensor)[1]   # Get rotation matrix
    rotated_molecule = new_coordinates.dot(np.linalg.inv(rotation))  # Rotate molecule to the new axis
    return 

# If nothing is specified, it will return the principal moments of inertia and the type of molecule
# If interested in other variable, write its name next to return

# Options: coordinates, atoms, atom_mass, mass_center, new_coordinates, tensor, diagonal,
#          I_i (i can be A, B or C), rotation and rotated_molecule
if __name__ == '__main__':
    molecule = input('Choose a .xyz file: \n')
    inertia(molecule)

