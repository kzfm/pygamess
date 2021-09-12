import pygamess
from pygamess import GamessError
from rdkit import Chem
import pytest


def test_ethane_sto3g():
    mol = Chem.MolFromMolFile("examples/ethane.mol", removeHs=False)
    g = pygamess.Gamess()
    g.debug = True
    try:
        r = g.run(mol)
    except GamessError:
        print(GamessError.value)

    assert pytest.approx(-78.30530748, 0.000000005) == r.mol.GetDoubleProp("total_energy")
