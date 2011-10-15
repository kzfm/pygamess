import sys
sys.path.append("../")
import pygamess

def test_gamess_ok():
    g = pygamess.Gamess(gamess_path="/usr/local/bin/gamess")
    assert g.gamess != None
