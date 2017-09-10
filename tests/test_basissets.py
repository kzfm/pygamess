from nose.tools import *
import pygamess
from pygamess import GamessError
from rdkit import Chem


def test_ethane_sto3g():
    mol = Chem.MolFromMolFile("examples/ethane.mol", removeHs=False)
    g = pygamess.Gamess(gamess_path="/usr/local/gamess")
    g.debug = True
    try:
        newmol = g.run(mol)
    except GamessError, gerr:
        print gerr.value

    eq_(newmol.GetDoubleProp("total_energy"), -78.30530748)
