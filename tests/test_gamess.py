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
    correct_input = ' $contrl scftyp=rhf runtyp=energy maxit=200 $end\n $basis gbasis=sto ngauss=3 $end\n $scf dirscf=.f. $end\n $system mwords=300 memddi=0 $end\n $DATA\n6324\nC1\nC      6.0     -0.7560000000    0.0000000000    0.0000000000 \nC      6.0      0.7560000000    0.0000000000    0.0000000000 \nH      1.0     -1.1404000000    0.6586000000    0.7845000000 \nH      1.0     -1.1404000000    0.3501000000   -0.9626000000 \nH      1.0     -1.1405000000   -1.0087000000    0.1781000000 \nH      1.0      1.1404000000   -0.3501000000    0.9626000000 \nH      1.0      1.1405000000    1.0087000000   -0.1781000000 \nH      1.0      1.1404000000   -0.6586000000   -0.7845000000 \n $END\n'
    assert g.input(mol) == correct_input

def test_parse_error():
    mol = Chem.MolFromMolFile("examples/ethane.mol", removeHs=False)
    g = pygamess.Gamess()
    g.charge(1)
    with pytest.raises(Exception) as e:
        r = g.run(mol)
    assert str(e.value).find("BUT YOU SELECTED MULTIPLICITY MULT=  1") > 0
