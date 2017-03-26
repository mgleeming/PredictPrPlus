#!/usr/bin/env python

from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
from Cython.Compiler.Options import directive_defaults

setup(
    name = 'predictPrPlus',
    version = '1.0.0',
    author = 'Michael Leeming',
    author_email = 'm.leeming@student.unimelb.edu.au',
    url = 'https://github.com/mgleeming/PredictPrPlus',
    license = 'LICENSE.txt',
    description = 'predictPrPlus: Predict protonation sites of multiply charged protein ions',

    packages = find_packages(),
    include_package_data = True,

    ext_modules = cythonize('./predictPrPlus/predictor/CSP.pyx'),

    entry_points = {
        'gui_scripts' : ['predictPrPlus = predictPrPlus.predictPrPlus:main']
    },

    install_requires = [
        'numpy',
        'pyqtgraph'
    ],
    
    dependency_links = ['http://pyqt.sourceforge.net/Docs/PyQt4/installation.html'],
    zip_safe = False
)
