#!/usr/bin/env python
#
# This file is subject to the terms and conditions defined in file
# 'LICENSE.txt', which is part of this source code package.

from setuptools import setup, find_packages

setup(
    name='statutils',
    version='0.2',
    description='Useful statistical routines and data plotting functions',
    author='Eric Suh',
    author_email='contact@ericsuh.com',
    packages=['statutils'],
    install_requires = [
        'scipy >= 0.10.1',
        'numpy >= 1.6.2',
        'matplotlib >= 1.2.0',
        'statsmodels >= 0.5.0',
    ],
    url='http://github.com/ericsuh/statutils',
    download_url='https://github.com/ericsuh/statutils/zipball/master',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Operating System :: OS Independent',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
    ],

    long_description="""\
Statutils
--------------

Collection of useful statistical routines and data plotting functions.
"""
)
