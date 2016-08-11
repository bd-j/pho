#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import re
import glob

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup

setup(
    name="pho",
    url="https://github.com/bd-j/pho",
    version='0.1.2',
    author="Ben Johnson",
    author_email="benjamin.johnson@cfa.harvard.edu",
    packages=["pho"],
    license="LICENSE",
    description="Large aperture photometry",
    long_description=open("README.md").read(),
    package_data={"": ["README.md", "LICENSE"]},
    scripts=[],
    include_package_data=True,
    install_requires=[],
)
