# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
import sys, os

here = os.path.abspath(os.path.dirname(__file__))
README = open(os.path.join(here, 'README.rst')).read()
NEWS = open(os.path.join(here, 'NEWS.txt')).read()

version = '0.3.3b1'

install_requires = [ ]

setup(name = 'pyXSteam',
      version=version,
      description = 'pyXSteam is a port of the Matlab/Excel Package XSteam by Magnus Holmgren, www.x-eng.com to Python',
      long_description = README + '\n\n' + NEWS,
      classifiers = ['Development Status :: 3 - Alpha',
                     'Programming Language :: Python',
                     'Topic :: Scientific/Engineering :: Physics'],
      keywords='steam water ice XSteam',
      author = 'drunsinn',
      author_email = 'dr.unsinn@googlemail.com',
      url = 'https://github.com/drunsinn/pyXSteam',
      license = 'LICENSE.txt',
      packages=find_packages(exclude=['tests*']),
      package_dir = {'pyXSteam': 'pyXSteam'},
      include_package_data=True,
      zip_safe = True,
      install_requires=install_requires,
      platforms = ('Any',),
      scripts = ['bin/pyXSteamDemo.py'],
      test_suite = 'pyXSteamTest.suite',
      tests_require = ['numpy >=1.6.2', ],
      use_2to3 = True,
)
