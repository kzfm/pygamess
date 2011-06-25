from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='pygamess',
      version=version,
      description="GAMESS wrapper for Python",
      long_description="""\
GAMESS wrapper for Python""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='chemistry',
      author='Ohkawa Kazufumi',
      author_email='kerolinq@gmail.com',
      url='https://github.com/kzfm/pygamess',
      license='MIT',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests', 'docs']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
