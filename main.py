'''
Python code to calculate SASA from given PDB and then calculate surface properties
using HPS-Urry or user provided parameters 
Author: Roshan M Regy
Email ID: roshanm.regy@tamu.edu 
''' 
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import argparse 
import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np

aa_param='''#AA  Charge  Lambda 
ALA       0.00      0.602942  
ARG       1.00      0.558824
ASN       0.00      0.588236
ASP       -1.00     0.294119
CYS       0.00      0.64706
GLN       0.00      0.558824  
GLU       -1.00     0.0
GLY       0.00      0.57353
HIS       0.00      0.764707  
ILE       0.00      0.705883
LEU       0.00      0.720589
LYS       1.00      0.382354
MET       0.00      0.676471 
PHE       0.00      0.82353 
PRO       0.00      0.758824
SER       0.00      0.588236 
THR       0.00      0.588236 
TRP       0.00      1.0 
TYR       0.00      0.897059 
VAL       0.00      0.664707'''

aa={}
for i in aa_param.split('\n'):
	if i[0]!='#':
		name=i.rsplit()[0]
		other=np.array(i.rsplit()[1:],dtype=float)
		aa[name]=other

# from Get Area http://curie.utmb.edu/area_man.html
random_coil_values = {
        "ALA":64.9,
        "ARG":195.5,
        "ASN":114.3,
        "ASP":113.0,
        "CYS":102.3,
        "GLN":143.7,
        "GLU":141.2,
        "HIS":154.6,
        "ILE":147.3,
        "GLY":87.2,
        "LEU":146.2,
        "LYS":164.5,
        "MET":158.3,
        "PHE":180.1,
        "PRO":105.2,
        "SER":77.4,
        "THR":106.2,
        "TRP":224.6,
        "TYR":193.1,
        "VAL":122.3
        }

# Use Lambert Azimuthal equal area projection for 
# mapping 3D coordinates to a 2D surface 
def lambert_azimuthal(coord):
    for i in range(coord.shape[1]):
        coord[:,i] = coord[:,i]/(np.sqrt(3)*np.max(np.absolute(coord[:,i])))
    zinv = np.sqrt(2/(1-coord[:,-1]))
    X = zinv*coord[:,0]
    Y = zinv*coord[:,1]
    return X,Y

# Bin 2D coordinates with weights from residue level properties
def bin2grid(X,Y,bins,weights):
    print ("X max: ",np.max(X)," X min: ",np.min(X))
    print ("Y max: ",np.max(Y)," Y min: ",np.min(Y))
    hist,xbinedge,ybinedge = np.histogram2d(X,Y,[bins,bins],weights=weights)
    return hist,xbinedge,ybinedge

# Get exposed residues from the sasa calculation and random coil values
# from GetArea  
def getexposed(m,c,chain):
    df = pd.DataFrame()
    nr = len(chain)
    sasa = np.zeros(nr)
    normsasa = np.zeros(nr)
    expocoord = []
    name = []
    weights = []
    for r,residue in enumerate(chain):
        name.append(residue.get_resname())
        print (name[-1])
        sasa[r] = residue.sasa
        normsasa[r] = sasa[r]/random_coil_values[name[-1]]
        if normsasa[r] > 0.5:
            print ("%s, %s, Exposed..."%(name[-1],normsasa[r]))
            expocoord.append(residue['CA'].get_coord())
            weights.append(aa[name[-1]])
    df['Name'] = name
    df['SASA'] = sasa
    df['NormSASA'] = normsasa
    df.to_csv('model%s_chain%s.csv'%(m,c))
    expocoord = np.array(expocoord)
    weights = np.array(weights)
    return expocoord, weights

# Main execution starts here
p = PDBParser(QUIET=1)
struct = p.get_structure('protein','sample.pdb')
sr = ShrakeRupley()
sr.compute(struct, level="R")
bins=10
for m,model in enumerate(struct):
    for c,chain in enumerate(model):
        expocoord, weights = getexposed(m,c,chain)
        X,Y = lambert_azimuthal(expocoord)
        for w in range(weights.shape[1]):
            fig,ax = plt.subplots(1,1,dpi=300,figsize=[3,3])
            hist,xbinedge,ybinedge = bin2grid(X,Y,bins,weights[:,w])
            im = ax.imshow(hist,cmap = "seismic")
            cb = plt.colorbar(im)
            plt.tight_layout()
            plt.savefig('imshow_%s.png'%w,dpi=300,transparent=True)
            plt.show()

