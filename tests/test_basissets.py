from nose.tools import *
import pygamess
import openbabel as ob

def test_ethane_sto3g():
    obc = ob.OBConversion()
    obc.SetInFormat("mol")
    mol = ob.OBMol()
    obc.ReadFile(mol, "../examples/ethane.mol")
    g = pygamess.Gamess(gamess_path="/usr/local/bin/rungms")
    newmol = g.run(mol)

    eq_(newmol.GetEnergy(),-78.30530748)

