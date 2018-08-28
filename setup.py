#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyXSteam",
    version="0.4.0",
    author="drunsinn",
    author_email="dr.unsinn@googlemail.com",
    keywords="steam water ice XSteam",
    description="pyXSteam is a port of the Matlab/Excel Package XSteam by Magnus Holmgren, www.x-eng.com to Python 3",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/drunsinn/pyXSteam",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Development Status :: 4 - Beta"
    ],
    packages=setuptools.find_packages(exclude=['docs', 'tests*']),
    install_requires=[],
    python_requires='>=3.6',
    include_package_data=True,
    zip_safe=True,
    scripts=['bin/pyXSteamDemo.py'],
    test_suite='pyXSteamTest.suite',
    tests_require=['numpy>=1.6.2', 'matplotlib>=2.2.3'],

)
