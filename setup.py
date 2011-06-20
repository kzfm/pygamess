from setuptools import setup, find_packages
import sys, os

version = '0.0'

setup(name='gamess',
      version=version,
      description="GAMESS wrapper for Python",
      long_description="""\
GAMESS wrapper for Python""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='chemistry',
      author='Ohkawa Kazufumi',
      author_email='kerolinq@gmail.com',
      url='http://blog.kzfmix.com',
      license='',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
