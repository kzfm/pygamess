from setuptools import setup, find_packages
import sys, os

version = '0.2.1'

setup(name='pygamess',
      version=version,
      description="GAMESS wrapper for Python",
      long_description="""\
`pygamess` is a GAMESS wrapper for Python

Requirements
------------
* Python 2.6 or later (not support 3.x)
* openbabel 2.3 
* GAMESS (and rungms script)

Features
--------
* nothing

Setup
-----
::

   $ easy_install pygamess

Basic Usage
-----------

single point calculation with pybel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

    >>> g.run_type('optimize')
    >>> optimized_mol = g.run(mol)
    >>> optimized_mol.GetEnergy()
    -78.306179642000004

changing basis-sets
~~~~~~~~~~~~~~~~~~~

basis_type method::

    >>> g.run_type('energy')
    >>> g.basis_type('sto3g')
    {'gbasis': 'sto', 'ngauss': '3'}
    >>> mol_sto3g = g.run(mol)
    >>> mol_sto3g.GetEnergy()
    -78.305307479999996
    >>> g.basis_type('631g')
    {'gbasis': 'N31', 'ndfunc': '1', 'ngauss': '6'}
    >>> mol_631g = g.run(mol)
    >>> mol_631g.GetEnergy()
    -79.228127109699997
    >>> g.basis_type('631gdp')
    {'gbasis': 'N31', 'ndfunc': '1', 'npfunc': '1', 'ngauss': '6'}
    >>> mol_631gdp = g.run(mol)
    >>> mol_631gdp.GetEnergy()
    -79.237634701499999

or edit property of Gamess instance::

    >>> g.basis = {'gbasis': 'sto', 'ngauss': '3'}
    >>> mol_sto3g = g.run(mol)
    >>> mol_sto3g.GetEnergy()
    -78.305307479999996

print GAMESS INPUT
~~~~~~~~~~~~~~~~~~

use gamess_input method

History
-------

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
