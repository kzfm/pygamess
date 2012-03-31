from setuptools import setup, find_packages
import sys, os

version = '0.3.0'

setup(name='pygamess',
      version=version,
      description="GAMESS wrapper for Python",
      long_description="""\
`pygamess` is a GAMESS wrapper for Python

Requirements
------------
* Python 2.6 or later (not support 3.x)
* openbabel 2.3 
* GAMESS

Features
--------
* nothing

Setup
-----
::

    $ easy_install pygamess

set GAMESS_PATH environment in your .bashrc::

    $ export GAMESS_HOME=/usr/local/gamess

Basic Usage
-----------

single point calculation with pybel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    >>> from pygamess import Gamess
    >>> import pybel
    >>> g = Gamess()
    >>> mol = pybel.readstring('smi','C=C')
    >>> mol.make3D()
    >>> optimized_mol = g.run(mol)
    >>> optimized_mol.energy
    -77.0722784993

single point calculation with openbabel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    >>> import pygamess
    >>> g = pygamess.Gamess()
    >>> import openbabel as ob
    >>> obc = ob.OBConversion()
    >>> obc.SetInFormat("mol")
    True
    >>> mol = ob.OBMol()
    >>> obc.ReadFile(mol, "examples/ethane.mol")
    True
    >>> try:
    ...     newmol = g.run(mol)
    ... except GamessError, gerr:
    ...     print gerr.value
    ... 
    >>> newmol.GetEnergy()
    -78.305307479999996

optimize calculation with pybel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set run_type::

    >>> from pygamess import Gamess
    >>> import pybel
    >>> g = Gamess()
    >>> g.run_type('optimize')
    >>> mol = pybel.readstring('smi','C')
    >>> mol.make3D()
    >>> optimized_mol = g.run(mol)
    >>> optimized_mol.energy
    -37.0895866208

changing basis-sets
~~~~~~~~~~~~~~~~~~~

basis_type method::

    >>> import pybel
    >>> from pygamess import Gamess
    >>> g = Gamess()
    >>> mol = pybel.readstring('smi', 'O')
    >>> mol.make3D()
    >>> g.run(mol).energy
    -74.96450135
    >>> g.run_type('optimize')
    >>> g.run(mol).energy
    -74.9659012146
    >>> g.basis_set('3-21G')
    {'gbasis': 'N21', 'ngauss': '3'}
    >>> g.run(mol).energy
    -75.585959758
    >>> g.basis_set('6-31G')
    {'gbasis': 'N31', 'ngauss': '6'}
    >>> g.run(mol).energy
    -75.9853591564
    >>> g.basis_set('6-311G')
    {'gbasis': 'N311', 'ngauss': '6'}
    >>> g.run(mol).energy
    -76.0109546389
    >>> g.basis_set('6-31G*')
    {'gbasis': 'N31', 'ndfunc': '1', 'ngauss': '6'}
    >>> g.run(mol).energy
    -76.0107465155
    >>> g.basis_set('6-31G**')
    {'gbasis': 'N31', 'ndfunc': '1', 'npfunc': '1', 'ngauss': '6'}
    >>> g.run(mol).energy
    -76.0236150193

or edit property of Gamess instance::

    >>> g.basis = {'gbasis': 'sto', 'ngauss': '3'}
    >>> mol_sto3g = g.run(mol)
    >>> mol_sto3g.GetEnergy()
    -78.305307479999996

print GAMESS INPUT
~~~~~~~~~~~~~~~~~~

use input method::

    >>> g.input(mol)


History
-------

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

""",
      classifiers=[
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Development Status :: 2 - Pre-Alpha',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Programming Language :: Python'
        ], 
      keywords='chemistry',
      author='Ohkawa Kazufumi',
      author_email='kerolinq@gmail.com',
      url='https://github.com/kzfm/pygamess',
      license='MIT',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests', 'docs']),
      include_package_data=True,
      zip_safe=False,
      requires=['openbabel']
      )
