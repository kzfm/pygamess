import pygamess
from rdkit import Chem
import pytest


def test_gamess_ok():
    g = pygamess.Gamess()
    assert type(g.rungms) != None


def test_ob_ok():
    mol = Chem.MolFromMolFile("examples/ethane.mol")
    assert type(mol) == Chem.rdchem.Mol


def test_gamess_input():
    mol = Chem.MolFromMolFile("examples/ethane.mol", removeHs=False)
    g = pygamess.Gamess()
    correct_input = ' $contrl scftyp=rhf runtyp=energy $end\n $basis gbasis=sto ngauss=3 $end\n $system mwords=30 $end\n $DATA\n6324\nC1\nC      6.0     -0.7560000000    0.0000000000    0.0000000000 \nC      6.0      0.7560000000    0.0000000000    0.0000000000 \nH      1.0     -1.1404000000    0.6586000000    0.7845000000 \nH      1.0     -1.1404000000    0.3501000000   -0.9626000000 \nH      1.0     -1.1405000000   -1.0087000000    0.1781000000 \nH      1.0      1.1404000000   -0.3501000000    0.9626000000 \nH      1.0      1.1405000000    1.0087000000   -0.1781000000 \nH      1.0      1.1404000000   -0.6586000000   -0.7845000000 \n $END\n'
    assert g.input(mol) == correct_input
