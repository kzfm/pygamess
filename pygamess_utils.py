from rdkit import Chem
from rdkit.Chem import AllChem

def rdkit_optimize(molobj):
    if type(molobj) is str:
        mol = Chem.MolFromSmiles(molobj)
        mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol,maxIters=200)
    return mol