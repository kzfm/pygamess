from rdkit import Chem
from rdkit.Chem import AllChem
from .gamout_parser import GamessOut

def rdkit_optimize(molobj):
    if type(molobj) is str:
        mol = Chem.MolFromSmiles(molobj)
        mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol,maxIters=200)
    return mol


def sdf2gamout(sdffile):
    suppl = Chem.SDMolSupplier(sdffile, removeHs=False)
    rs = []
    for m in suppl:
        r = GamessOut()
        r.total_energy = m.GetDoubleProp("total_energy")
        r.HOMO = m.GetDoubleProp("HOMO")
        r.LUMO = m.GetDoubleProp("LUMO")
        r.nHOMO = m.GetDoubleProp("nHOMO")
        r.nLUMO = m.GetDoubleProp("nLUMO")
        r.dipole_moment = [
            m.GetDoubleProp("dx"),
            m.GetDoubleProp("dy"),
            m.GetDoubleProp("dz"),
            m.GetDoubleProp("dipole_moment")]
        r.orbital_energies = eval(m.GetProp("orbital_energies"))
        
        if m.HasProp("uv_spectra") > 0:
            r.uv_spectra = eval(m.GetProp("uv_spectra"))
        if m.HasProp("isotropic_shielding") > 0:
            r.isotropic_shielding = eval(m.GetProp("isotropic_shielding"))
        if m.HasProp("ir_spectra") > 0:
            r.ir_spectra = eval(m.GetProp("ir_spectra"))

        r.mulliken_charges = [a.GetDoubleProp("mulliken_charge") for a in m.GetAtoms()]
        r.lowdin_charges = [a.GetDoubleProp("lowdin_charge") for a in m.GetAtoms()]
        r.mulliken_populations = [a.GetDoubleProp("mulliken_population") for a in m.GetAtoms()]
        r.lowdin_populations = [a.GetDoubleProp("lowdin_population") for a in m.GetAtoms()]
        r.mol = m
        rs.append(r)

    return rs