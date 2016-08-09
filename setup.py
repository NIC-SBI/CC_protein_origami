"""
This is a basic setuptools file.
"""

import sys

from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

class PyTest(TestCommand):
    """
    Unit test wrapper for the PyTest, including coverage repport
    """
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['--color=yes', 'tests/'] #['--cov spectraplotpy tests/']
        self.test_suite = True
    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)

setup(
    name="ppmod",
    version="0.0.1",
    packages=find_packages('ppmod'),
    scripts=[],

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=[
        #"numpy", numpy is best installed via conda etc...
    ],

    package_data={
        # If any package contains *.txt or *.rst files, include them:
        '': 'ppmod',
        '': ['*.txt', '*.md'],
        '': 'examples',
        '': 'building-blocks',
        '': 'notebooks'
        # And include any *.msg files found in the 'hello' package, too:
        #'hello': ['*.msg'],
    },

    # Project uses pytest for the tests

    tests_require=[
        'pytest',
        'pytest-cov',
        'mock'
    ],
    
    cmdclass={
        'test': PyTest
    },

    # metadata for upload to PyPI
    author="Ajasja Ljubetic",
    author_email="ajasja.ljubetic@gmail.com",
    description="A package to build models of protein polyhedra",
    license="None",
    keywords="protein polyhedra",

    # project home page, if any
#    url="https://github.com/SpectralGroup/spectraplotpy"   

    # could also include long_description, download_url, classifiers, etc.
)