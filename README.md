With this program, it is possible to calculate the inertia tensor and the principal moments of inertia of a molecule (given in XYZ format). It also classifies the molecule as a prolate/oblate, spherical top or asymmetrical top.

First, we have a function that reads a .xyz file and returns a list of atoms and a N ×3 array with the coordinates of the atoms. The following function contains a dictionary with atoms up to Cl and transforms the atom list into a one-dimensional array with the respective atom masses. center of mass,as its name indicates, computes the center of mass to then subtract it from the original coordinates.


The next block is the most important one. I, generates the inertia tensor by taking as arguments the translated geometry and the masses. The Kronecker delta is programmed as a two-case if block to differentiate the diagonal and offdiagonal elements. diag takes the inertia tensor and returns several variables:

• diag(I(...))[0] = D = diagonal matrix with the eigenvalues sorted.

• diag(I(...))[1] = R = rotation matrix made of normalized eigenvectors
in columns.

• diag(I(...))[2/3/4] = Ia/Ib/Ic = principal moments of inertia.


The top function recognises the type of molecule according to the principal moments of inertia.


Lastly, inertia, groups all the previous functions and adds the rotated molecule variable, which contains the coordinates of the molecule oriented in the principalaxis direction. By default, this function only returns the type of molecule and the principal moments of inertia in g × cm2 units, however, it can be used to get any of the previously defined variables by adding the keyword next to return.

The program was tested with four molecules:

• Methane (CH4): spherical top

• Chloroform (CH3Cl): prolate symmetrical top

• Benzene (C6H6): oblate symmetrical top

• Water (H2O): asymmetrical top

In the case of CH4, the output is the following:
The molecule is a spherical top
I_A = 5.105345287105582e-40 gcm^2,
I_B = 5.105343554617158e-40 gcm^2,
I_C = 5.105343029325426e-40gcm^2
And according to literature [1], the principal inertia moment of methane is 5.332 × 10−40, so the result can be labeled as accurate.


For CH3Cl:
The molecule is a prolate symmetrical top
I_A = 6.186620168262775e-39 gcm^2,
I_B = 6.186620061792865e-39 gcm^2,
I_C = 5.105343746642018e-40gcm^2
The highest inertia moments’ experimental value (IAandIB) is 6.313×10−39 and for the other one (IC): 5.385 × 10−40 [1]. So again, they are close and can be treated as correct.


Now with benzene:
The molecule is an oblate symmetrical top
I_A = 2.9765008234292364e-38 gcm^2,
I_B = 1.4882511393384954e-38 gcm^2,
I_C = 1.4882496840907405e-38gcm^2
The highest experimental rotational constant (IA) is: 2.953×10−38, and the other two (IBandIC): 1.476 × 10−38 [1]. In this case, the inertia moments are higher due to the higher mass and size of the molecule.


Finally, for water:
The molecule is an asymetric top
I_A = 2.9674124065143283e-40 gcm^2,
I_B = 2.05553148718378e-40 gcm^2,
I_C = 9.118809193305478e-41gcm^2
The experimental values for water are: IA = 3.015 × 10−40, IB = 1.929 × 10−40 and IC = 1.004 × 10−40 [1].


With all this information, we can conclude that the program calculates correctly these principal moments of inertia, and the small deviations from the
experimental value can be rationalized by the fact that the geometries are not optimized.

[1] G. Herzberg. Electronic spectra and electronic structure of polyatomic molecules.New York. Van Nostrand, 1966
