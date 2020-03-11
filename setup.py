#!/usr/bin/env python
## -*- encoding: utf-8 -*-

import os
import sys
from setuptools import setup
from codecs import open # To open the README file with proper encoding
from setuptools.command.test import test as TestCommand # for tests


# Get information from separate files (README, VERSION)
def readfile(filename):
    with open(filename,  encoding='utf-8') as f:
        return f.read()

# For the tests
class SageTest(TestCommand):
    def run_tests(self):
        errno = os.system("sage -t --force-lib pycontrolledreduction")
        if errno != 0:
            sys.exit(1)






setup(
    name="crystalline_obstruction",
    author="Edgar Costa, Emre Sertoz",
    author_email="edgarc@mit.edu",
    url="https://github.com/edgarcosta/crystalline_obstruction",
    license="GNU General Public License, version 2 or 3",
    description="Wrapper for computing an approximation to the crystalline obstruction map ",
    long_description = readfile("README.rst"), # get the long description from the README
    version = readfile("VERSION"), # the VERSION file is shared with the documentation
    classifiers=[
      # How mature is this project? Common values are
      #   3 - Alpha
      #   4 - Beta
      #   5 - Production/Stable
      'Development Status :: 4 - Beta',
      'Intended Audience :: Science/Research',
      'Topic :: Scientific/Engineering :: Mathematics',
      'License :: OSI Approved :: GNU General Public License v2 or v3',
      'Programming Language :: Python :: 3.7',
    ], # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords = "sagemath crystalline obstruction",
    setup_requires=["sagemath"], # currently useless, see https://www.python.org/dev/peps/pep-0518/
    install_requires=["pycontrolledreduction", "sagemath", "sphinx"],
    packages=["crystalline_obstruction"],
    include_package_data = True,
    cmdclass = {'test': SageTest}, # adding a special setup command for tests
)
