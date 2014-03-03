# -*- coding: utf-8 -*-
from setuptools import setup


setup(
    name = 'pyXSteam',
    version = '0.2',
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
    url = '',
    license = 'LICENSE.txt',
    # requires = ['numpy >=1.6.2', ],
    test_suite = 'test.suite',
    tests_require = ['numpy >=1.6.2', ]
)

# sudo pip install .
# python setup.py bdist_egg
# python setup.py sdist
# python setup.py test
