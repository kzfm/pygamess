from nose.tools import *
import pygamess
import openbabel as ob

def test_gamess_ok():
    g = pygamess.Gamess(gamess_path="/usr/local/bin/rungms")
    assert g.gamess != None

def test_ob_ok():
    obc = ob.OBConversion()
    obc.SetInFormat("mol")
    mol = ob.OBMol()
    r = obc.ReadFile(mol, "examples/ethane.mol")
    assert r == 1

def test_gamess_input():
    obc = ob.OBConversion()
    obc.SetInFormat("mol")
    mol = ob.OBMol()
    obc.ReadFile(mol, "examples/ethane.mol")
    g = pygamess.Gamess(gamess_path="/usr/local/bin/rungms")
    print g.gamess_input(mol)

    correct_input = """ $contrl runtyp=energy scftyp=rhf mult=1  $end
 $basis gbasis=sto ngauss=3 $end
 $system mwords=30  $end
 $DATA
6324
C1
C      6.0     -0.7560000000    0.0000000000    0.0000000000 
C      6.0      0.7560000000    0.0000000000    0.0000000000 
H      1.0     -1.1404000000    0.6586000000    0.7845000000 
H      1.0     -1.1404000000    0.3501000000   -0.9626000000 
H      1.0     -1.1405000000   -1.0087000000    0.1781000000 
H      1.0      1.1404000000   -0.3501000000    0.9626000000 
H      1.0      1.1405000000    1.0087000000   -0.1781000000 
H      1.0      1.1404000000   -0.6586000000   -0.7845000000 
 $END


"""
    eq_(g.gamess_input(mol), correct_input)
