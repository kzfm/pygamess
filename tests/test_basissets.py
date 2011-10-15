import sys

sys.path.append("../")

import pygamess
import openbabel as ob
g = pygamess.Gamess()
obc = ob.OBConversion()
obc.SetInFormat("mol")
mol = ob.OBMol()
obc.ReadFile(mol, "../examples/ethane.mol")

def test_ethane_sto3g():
    newmol = g.run(mol)
    assert newmol.GetEnergy() == -78.305307479999996
