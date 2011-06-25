==========
 pygamess
==========

Basic Usage
===========

import module and read molecule from file::

    >>> import pygamess
    >>> import openbabel as ob
    >>> g = pygamess.Gamess()
    >>> obc = ob.OBConversion()
    >>> obc.SetInFormat("mol")
    True
    >>> mol = ob.OBMol()
    >>> obc.ReadFile(mol, "examples/ethane.mol")
    True

run GAMESS::

    >>> try:
    ...     newmol = g.run(mol)
    ... except GamessError, gerr:
    ...     print gerr.value
    ... 

get energy::

    >>> newmol.GetEnergy()
    -78.305307479999996
