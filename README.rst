==========
 pygamess
==========

`pygamess` is a GAMESS wrapper for Python

Requirements
------------
* Python 2.6 or later (not support 3.x)
* RDKit
* GAMESS

Features
--------
* nothing

Setup
-----
::

    $ pip install pygamess

set GAMESS_HOME environment in your .bashrc::

    $ export GAMESS_HOME=/usr/local/gamess

Basic Usage
-----------

single point calculation with RDKit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    >>> from pygamess import Gamess
    >>> from rdkit import Chem
    >>> from rdkit.Chem import AllChem
    >>> m = Chem.MolFromSmiles("CC")
    >>> m = Chem.AddHs(m)
    >>> AllChem.EmbedMolecule(m)
    0
    >>> AllChem.UFFOptimizeMolecule(m,maxIters=200)
    0
    >>> g = Gamess()
    >>> nm = g.run(m)
    >>> nm.GetProp("total_energy")
    '-78.302511990200003'

optimize calculation with RDKit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Change run_type::

    >>> from pygamess import Gamess
    >>> from rdkit import Chem
    >>> from rdkit.Chem import AllChem
    >>> m = Chem.MolFromSmiles("CC")
    >>> m = Chem.AddHs(m)
    >>> AllChem.EmbedMolecule(m)
    0
    >>> AllChem.UFFOptimizeMolecule(m,maxIters=200)
    0
    >>> g = Gamess()
    >>> g.run_type('optimize')
    >>> optimized_mol = g.run(mol)
    >>> optimized_mol.GetProp("total_energy")
    '-78.306179626599999'

change basis-sets
~~~~~~~~~~~~~~~~~

Change basis set::

    >>> from pygamess import Gamess
    >>> from rdkit import Chem
    >>> from rdkit.Chem import AllChem
    >>> m = Chem.MolFromSmiles("CC")
    >>> m = Chem.AddHs(m)
    >>> AllChem.EmbedMolecule(m)
    0
    >>> AllChem.UFFOptimizeMolecule(m,maxIters=200)
    0
    >>> g = Gamess()
    >>> g.run_type = "optimize"
    >>> g.run(m).GetProp("total_energy")
    '-78.302511907400003'
    >>> g.basis_set("3-21G")
    {'gbasis': 'N21', 'ngauss': '3'}
    >>> g.run(m).GetProp("total_energy")
    '-78.790442552000002'
    >>> g.basis_set("6-31G")
    {'gbasis': 'N31', 'ngauss': '6'}
    >>> g.run(m).GetProp("total_energy")
    '-79.194024566899998'
    >>> g.basis_set("6-31G*")
    {'gbasis': 'N31', 'ndfunc': '1', 'ngauss': '6'}
    >>> g.run(m).GetProp("total_energy")
    '-79.225521673399996'
    >>> g.basis_set("6-31G**")
    {'gbasis': 'N31', 'ndfunc': '1', 'npfunc': '1', 'ngauss': '6'}
    >>> g.run(m).GetProp("total_energy")
    '-79.235082450899995'
    
or edit property of Gamess instance::

    >>> g.basis = {'gbasis': 'sto', 'ngauss': '3'}
    >>> mol_sto3g = g.run(m)
    >>> mol_sto3g.GetProp("total_energy")
    '-78.302511907400003'

print GAMESS INPUT
~~~~~~~~~~~~~~~~~~

use input method::

    >>> g.input(mol)


History
-------

0.4.0 (2017-09-13)
~~~~~~~~~~~~~~~~~~

* Changed the backend library from openbabel to RDKit

0.3.0 (2012-03-31)
~~~~~~~~~~~~~~~~~~

* no more required rungms script and use internal rungms (default)
* added basis_set method(STO-3G,3-21G,6-31G,6-311G,6-31G*,6-31G**,AM1,PM3,MNDO)
* constructor can accept options
* bug fixed (spin multiplicity)

0.2.2 (2012-03-30)
~~~~~~~~~~~~~~~~~~

* added charge settings
* method name changed (gamess_input -> input)

0.2.1 (2012-03-23)
~~~~~~~~~~~~~~~~~~

* bug fixed (multiplicity setting for pybel) 
* bug fixed (print error when rungms exec failed)
* added document

0.2.0 (2012-03-06)
~~~~~~~~~~~~~~~~~~

* run method can accept OBMol and Pybel-Molecule object

0.1.2 (2011-09-23)
~~~~~~~~~~~~~~~~~~

* added CIS method (and optimization)

0.1.1 (2011-08-06)
~~~~~~~~~~~~~~~~~~

* updated document
* semiempical method (AM1, PM3, MNDO)
* added statpt option
* changed default error print (10 lines)

0.1 (2011-6-25)
~~~~~~~~~~~~~~~~~~
* first release
