#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 12:24:58 2022

@author: hagen
"""

import sys

required_verion = (3,)
if sys.version_info < required_verion:
    raise ValueError('SurfRadPy needs at least python {}! You are trying to install it under python {}'.format('.'.join(str(i) for i in required_verion), sys.version))

from setuptools import setup, find_packages

setup(
    name="lrt_wrap",
    version="0.1.1", #setting this caused a huge hadeache ... basically the script wasn't found when the version was set
    packages=find_packages(),
    author="Hagen Telg",
    author_email="hagen.telg@noaa.gov",
    description="...",
    license="MIT",
    url="https://github.com/hagne/lrt_wrap",
    # install_requires=['pandas', 'numpy', 'xarray'],
    # scripts=['scripts/scrape_hrrr', 'scripts/modelextractor',
    #          # 'scripts/hrrr_smoke2gml'
    #          ],
    # entry_points = {'console_scripts': ['qcrad2ncei=SurfRadPy.NCEI:qcrad2ncei'],},
    # package_data={'': ['*.cdl']},
    # include_package_data=True,
    # zip_safe=False
)