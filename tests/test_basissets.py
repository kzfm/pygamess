import pygamess
import openbabel as ob

def test_ethane_sto3g():
    obc = ob.OBConversion()
    obc.SetInFormat("mol")
    mol = ob.OBMol()
    obc.ReadFile(mol, "examples/ethane.mol")
    g = pygamess.Gamess(gamess_path="/usr/local/bin/rungms")
    print g.gamess_input(mol)
    newmol = g.run(mol)

    assert newmol.GetEnergy() == -78.30530748

if __name__ == '__main__':
    test_ethane_sto3g()
