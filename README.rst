==========
 pygamess
==========

`pygamess` is a GAMESS wrapper for Python


Requirements
------------
* Python 3.7 or later (pygamess <0.5 supports only Python2)
* RDKit >= 2020.03.5
* GAMESS
* ruamel.YAML

Setup
-----
::

    $ pip install pygamess

set GAMESS_HOME environment in your .bashrc or .zshrc::

    $ export GAMESS_HOME=/usr/local/gamess

Test
-----
::

    $ pytest

Basic Usage
-----------

Single point calculation
~~~~~~~~~~~~~~~~~~~~~~~~
::

    >>> from pygamess import Gamess
    >>> from rdkit import Chem
    >>> from pygamess_utils import rdkit_optimize
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess()
    >>> r = g.run(m)
    >>> r.total_energy
    -152.127991054

Or use rdkit directly::

    >>> from pygamess import Gamess
    >>> from rdkit import Chem
    >>> from rdkit.Chem import AllChem
    >>> m = Chem.MolFromSmiles("CCO")
    >>> m = Chem.AddHs(m)
    >>> AllChem.EmbedMolecule(m)
    0
    >>> AllChem.UFFOptimizeMolecule(m,maxIters=200)
    0
    >>> g = Gamess()
    >>> r = g.run(m)
    >>> r.total_energy
    -152.1279910526

The GamessOut object(r) contains the results of the GAMESS calculation 
and the RDKit Chem object with the calculation results::

    >>> r.total_energy
    -152.1279910526
    >>> r.HOMO
    -0.3453
    >>> r.nHOMO
    -0.3978
    >>> r.LUMO
    0.5594
    >>> r.nLUMO
    0.6127
    >>> r.dipole_moment
    [0.681619, -0.605188, 1.146253, 1.464497]
    >>> r.orbital_energies
    [-20.2521, -11.0932, -11.0402, -1.286, -0.9614, -0.7909, -0.6269, -0.5716, -0.5347, -0.4976, -0.4705, -0.3978, -0.3453, 0.5594, 0.6127, 0.6639, 0.69, 0.7002, 0.7388, 0.7549, 0.7852]
    >>> r.mulliken_charges
    [-0.171226, 0.024351, -0.298162, 0.049615, 0.055831, 0.061924, 0.042714, 0.061112, 0.173842]
    >>> r.lowdin_charges
    [-0.089262, 0.062464, -0.212689, 0.022024, 0.02752, 0.031644, 0.010235, 0.026611, 0.121453]

The RDKit Chem object has the same information to store these data into SDF::

    >>> r.mol
    <rdkit.Chem.rdchem.Mol object at 0x7ffd48217f30>
    >>> r.mol.GetProp("total_energy")
    '-152.12799105260001'
    >>> r.mol.GetProp("HOMO")
    '-0.3453'
    >>> r.mol.GetProp("nHOMO")
    '-0.39779999999999999'
    >>> r.mol.GetProp("LUMO")
    '0.55940000000000001'
    >>> r.mol.GetProp("nLUMO")
    '0.61270000000000002'
    >>> r.mol.GetProp("dipole_moment")
    '1.4644969999999999'
    >>> r.mol.GetProp("dx")
    '0.68161899999999997'
    >>> r.mol.GetProp("dy")
    '-0.60518799999999995'
    >>> r.mol.GetProp("dz")
    '1.146253'
    >>> r.mol.GetProp("orbital_energies")
    '-20.2521 -11.0932 -11.0402 -1.286 -0.9614 -0.7909 -0.6269 -0.5716 -0.5347 -0.4976 -0.4705 -0.3978 -0.3453 0.5594 0.6127 0.6639 0.69 0.7002 0.7388 0.7549 0.7852'
    >>> for a in r.mol.GetAtoms():
    ...   print("{}:\t{:.4f}\t{:.4f}".format(a.GetSymbol(), float(a.GetProp("mulliken_charge")), float(a.GetProp("lowdin_charge"))))
    ... 
    C:	-0.1712	-0.0893
    C:	0.0244	0.0625
    O:	-0.2982	-0.2127
    H:	0.0496	0.0220
    H:	0.0558	0.0275
    H:	0.0619	0.0316
    H:	0.0427	0.0102
    H:	0.0611	0.0266
    H:	0.1738	0.1215


Geometry optimization
~~~~~~~~~~~~~~~~~~~~~

Set the run_type as 'optimize'. This optimization process updates the coordinates of the molecule::

    >>> from pygamess import Gamess
    >>> from pygamess_utils import rdkit_optimize
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess()
    >>> g.run_type('optimize')
    >>> r = g.run(mol)
    >>> r.total_energy
    -152.1330661028
    >>> original_conf = m.GetConformer(0)
    >>> optimized_conf = r.mol.GetConformer(0)
    >>> for c in original_conf.GetPositions():
    ...   print(c)
    ... 
    [ 0.91206647 -0.11944851 -0.1294722 ]
    [-0.47153193  0.42043351  0.21118521]
    [-1.44831334 -0.21539324 -0.56715297]
    [ 0.9650486  -1.2050043   0.09903891]
    [1.67955732 0.41189183 0.47186063]
    [ 1.12654515  0.0378618  -1.20780162]
    [-0.67995345  0.27545635  1.29552482]
    [-0.49851483  1.50993204 -0.00177488]
    [-1.58490399 -1.11572948 -0.17140791]
    >>> for c in optimized_conf.GetPositions():
    ...   print(c)
    ... 
    [ 0.91442972 -0.13086468 -0.12174822]
    [-0.48373921  0.42850882  0.23169745]
    [-1.54145595 -0.17763397 -0.52424945]
    [ 0.97874385 -1.18768306  0.12185969]
    [1.67944907 0.39674369 0.43935649]
    [ 1.11452802 -0.00923776 -1.18177409]
    [-0.65534535  0.31918718  1.31002371]
    [-0.51450285  1.49751073  0.00304491]
    [-1.4921073  -1.13653094 -0.2782105 ]

Changing basis sets
~~~~~~~~~~~~~~~~~~~

Use basis_sets method::

    >>> from pygamess import Gamess
    >>> from pygamess_utils import rdkit_optimize
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess()
    >>> g.run_type = "optimize"
    >>> g.run(m).total_energy
    -152.127991054
    >>> g.basis_sets("3-21G")
    >>> g.run(m).total_energy
    -153.2170653562
    >>> g.basis_sets("6-31G")
    >>> g.run(m).total_energy
    -154.0054866151
    >>> g.basis_sets("6-31G*")
    >>> g.run(m).total_energy
    -154.0702703669
    >>> g.basis_sets("6-31G**")
    >>> g.run(m).total_energy
    -154.0843823698

Or edit the basis attribute directly::

    >>> g.options['basis'] = {'gbasis': 'sto', 'ngauss': '3'}
    >>> g.run(m).total_energy
    -152.127991054

DFT calculation
~~~~~~~~~~~~~~~

B3LYP/6-31G*::

    >>> from pygamess import Gamess
    >>> from pygamess_utils import rdkit_optimize
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess()
    >>> g.run_type("optimize")
    >>> g.basis_sets("6-31G*")
    >>> g.dft_type("B3LYP")
    >>> g.run(m).total_energy
    -154.9387962055


M062X/6-31G**::

    >>> from pygamess import Gamess
    >>> from pygamess_utils import rdkit_optimize
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess()
    >>> g.run_type("optimize")
    >>> g.dft_type("M06-2X")
    >>> g.basis_sets("6-31G**")
    >>> g.run(m).total_energy
    -154.9636095207

PCM calculation
~~~~~~~~~~~~~~~

Pygamess currently only supports CPCM, but will support IEFPCM in the future::

    >>> from pygamess import Gamess
    >>> from pygamess_utils import rdkit_optimize
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess()
    >>> g.basis_sets("6-31G*")
    >>> g.pcm_type("water")
    >>> g.run_type("optimize")
    >>> r = g.run(m)
    >>> r.total_energy
    -154.0824604616
    >>> r.internal_energy
    -154.0748367584
    >>> r.delta_internal_energy
    0.0
    >>> r.electrostatic_interaction
    -0.0076237032
    >>> r.pierotti_cavitation_energy
    0.0
    >>> r.dispersion_free_energy
    0.0
    >>> r.repulsion_free_energy
    0.0
    >>> r.total_interacion
    -0.0076237032

Printing GAMESS INPUT
~~~~~~~~~~~~~~~~~~~~~

use input method::

    >>> from pygamess import Gamess
    >>> from pygamess_utils import rdkit_optimize
    >>> m = rdkit_optimize("CO")
    >>> g = Gamess()
    >>> print(g.input(m))
    $contrl scftyp=rhf runtyp=energy $end
    $basis gbasis=sto ngauss=3 $end
    $system mwords=100 $end
    $DATA
    6324
    C1
    C      6.0     -0.3577002260    0.0075902163   -0.0214817423 
    O      8.0      0.9087355734   -0.5349924519   -0.2611189822 
    H      1.0     -0.5468334701    0.0717914414    1.0721087268 
    H      1.0     -0.4337681128    1.0193437527   -0.4757947304 
    H      1.0     -1.1269974200   -0.6479305528   -0.4789564639 
    H      1.0      1.5565636556    0.0841975943    0.1652431920 
    $END

Debugging pygamess
~~~~~~~~~~~~~~~~~~

set PYGAMESS_DEBUG environment::

    $ export PYGAMESS_DEBUG=1

This won't remove the all files generated by the GAMESS executable, including the output files.

set logger level::

    >>> from pygamess import Gamess, logger
    >>> import logging
    >>> logger.setLevel(logging.DEBUG)
    >>> g = Gamess()
    DEBUG:pygamess.gamess:tmpdir: /var/folders/gm/4tcnnyqd09d2jt7p0dtvr28m0000gn/T/tmp889j9c7e

History
-------

0.6.0 (2021-05-02)
~~~~~~~~~~~~~~~~~~~~

* Support DFT calculation
* Support PCM calculation (C-PCM only)
* Improve the parser
* Support logger levels
* Change method name from "basis_set" to "basis_sets"

0.5.0 (2020-09-13)
~~~~~~~~~~~~~~~~~~~~

* Support Python3

0.4.1.1 (2017-09-16)
~~~~~~~~~~~~~~~~~~~~

* Update Readme

0.4.1 (2017-09-16)
~~~~~~~~~~~~~~~~~~

* Bug fix (coordinates problem)

0.4.0 (2017-09-13)
~~~~~~~~~~~~~~~~~~

* Change the backend library from openbabel to RDKit

0.3.0 (2012-03-31)
~~~~~~~~~~~~~~~~~~

* Use internal rungms (default)
* Add basis_set method(STO-3G,3-21G,6-31G,6-311G,6-31G*,6-31G**,AM1,PM3,MNDO)
* Constructor can accept options
* Bug fixed (spin multiplicity)

0.2.2 (2012-03-30)
~~~~~~~~~~~~~~~~~~

* Add charge settings
* Change Method name (gamess_input -> input)

0.2.1 (2012-03-23)
~~~~~~~~~~~~~~~~~~

* Bug fix (multiplicity setting for pybel) 
* Bug fix (print error when rungms exec failed)
* Add document

0.2.0 (2012-03-06)
~~~~~~~~~~~~~~~~~~

* Run method accepts OBMol and Pybel-Molecule object

0.1.2 (2011-09-23)
~~~~~~~~~~~~~~~~~~

* Add CIS method (and optimization)

0.1.1 (2011-08-06)
~~~~~~~~~~~~~~~~~~

* Update document
* Semiempical method (AM1, PM3, MNDO)
* Add statpt option
* Change default error print message (10 lines)

0.1 (2011-6-25)
~~~~~~~~~~~~~~~~~~
* First release
