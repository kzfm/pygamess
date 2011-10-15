import sys
sys.path.append("../")
import pygamess

def test_gamess_ok():
    g = pygamess.Gamess()
    assert g.gamess != None
