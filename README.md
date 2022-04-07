# surf2d
Plot surface properties of folded proteins onto a 2D map
The code uses Biopython, numpy and Matplotlib as prequisites.

The purpose of the code is to present surface features of folded proteins on a 2D map. The properties of the surface residues can be decided based on user preferences. Currently, I am using HPS-Urry parameters (hydrophobicity, charge) to calculate surface properties.

It follows the following steps,
1) use Bio.PDB to parse a PDB file
2) use the Shrake-Rupley algorithm to find exposed surface area of all residues in the PDB file
3) use random coil values from GetArea to normalize and classify residues as exposed or buried
4) Convert the C-alpha 3D coordinates of exposed residues to a 2D map using the Lambert Azimuthal equal area projection formula 
5) Bin these 2D coordinates into a 2D histogram and assign weights from the coarse grained model for each residue as it is binned.

Corrections to be made: 
1) I want to use Sanner's MSMS algorithm to find residue depth and use that measure to define the protein surface and compare with the normal rolling ball method.
2) Make the code also do these calculations on the atom level for atomistic simulations.
3) Try other methods of mapping the coordinates to a 2D surface.

The code is very rough at the moment and would change considerably.
