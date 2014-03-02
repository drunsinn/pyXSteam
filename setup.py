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
    author_email = 'pirkl.max@googlemail.com',
    classifiers = [ 'Development Status ::  Alpha',
                    'Programming Language :: Python',
                    ],
    platforms = ('Any',),
    packages = ['pyXSteam', 'test'],
    scripts = ['bin'],
    url = 'http://pypi.python.org/pypi/TowelStuff/',
    license = 'LICENSE.txt',
    requires = ['numpy(>=1.6.2)', ],
    test_suite = 'test.TestXSteam_MKS',
)

# sudo pip install .
# python setup.py bdist_egg
# python setup.py sdist
# python setup.py test
