import sys
sys.path.append("../")
import pygamess
import openbabel as ob

def test_ethane_sto3g():
    obc = ob.OBConversion()
    obc.SetInFormat("mol")
    mol = ob.OBMol()
    obc.ReadFile(mol, "../examples/ethane.mol")
    g = pygamess.Gamess(gamess_path="/usr/local/bin/rungms")
    print g.gamess_input(mol)
    try:
        newmol = g.run(mol)
    except GamessError, gerr:
        print gerr.value

    print  newmol.GetEnergy()
    assert 1 == 1

if __name__ == '__main__':
    test_ethane_sto3g()
