# -*- coding: utf-8 -*-
'''
Created on 03.02.2014

@author: max
'''
from setuptools import setup

setup(
    name = 'pyXSteam',
    version = '0.1',
    description = 'pyXSteam is a port of the Matlab/Excel Package XSteam by Magnus Holmgren, www.x-eng.com to Python',
    long_description = open('README.txt').read(),
    author = 'Max Pirkl',
    author_email = 'dr.unsinn@googlemail.com',
    classifiers = [ 'Development Status ::  Alpha',
                    'Programming Language :: Python',
                    ],
    platforms = ('Any',),
    packages = ['pyXSteam', 'pyXSteam.test'],
    url = 'http://pypi.python.org/pypi/TowelStuff/',
    license = 'LICENSE.txt',
    requires = ['numpy(>=1.6.2)', ],
    test_suite = 'pyXSteam.test.test_pyXSteam',
)

# sudo pip install .
# python setup.py bdist_egg
# python setup.py sdist
# python setup.py test
