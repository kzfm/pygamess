from nose.tools import *
import sys
sys.path.append("../")
import pygamess
import openbabel as ob

def test_ethane_sto3g():
    obc = ob.OBConversion()
    obc.SetInFormat("mol")
    mol = ob.OBMol()
    obc.ReadFile(mol, "../examples/ethane.mol")
    g = pygamess.Gamess(gamess_path="/usr/local/bin/gamess")
    try:
        newmol = g.run(mol)
    except GamessError, gerr:
        print gerr.value
    assert newmol.GetEnergy() == -78.305307479999996
