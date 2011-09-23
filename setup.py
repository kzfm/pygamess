from setuptools import setup, find_packages
import sys, os

version = '0.1.2'

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

History
-------

0.1.2 (2011-9-23)
~~~~~~~~~~~~~~~~~~
* added CIS method (and optimization)

0.1.1 (2011-8-6)
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
