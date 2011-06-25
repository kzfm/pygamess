from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='pygamess',
      version=version,
      description="GAMESS wrapper for Python",
      long_description="""\
`pygamess` is a GAMESS wrapper for Python

Requirements
------------
* Python 2.6 or later (not support 3.x)
* openbabel 2.3 

Features
--------
* nothing

Setup
-----
::

   $ easy_install pygamess

History
-------
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
