# setup.py script adapted from COMBAT sc_pipelines repo, originally by Adam Cribbs

import sys
import os
import re
import setuptools
from setuptools import setup, find_packages, Extension

from distutils.version import LooseVersion
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    print("Version detected:", LooseVersion(setuptools.__version__))
    raise ImportError(
        "sc_pipelines requires setuptools 1.1 higher")

########################################################################
########################################################################
IS_OSX = sys.platform == 'darwin'

########################################################################
########################################################################
# collect version
print(sys.path.insert(0, "scpipelines"))
import version

version = version.__version__

###############################################################
###############################################################
# Define dependencies
#
major, minor1, minor2, s, tmp = sys.version_info

if major < 3:
    raise SystemExit("""Requires Python 3 or later.""")

sc_pipelines_packages = find_packages()
sc_pipelines_package_dirs = {'pipelines': 'pipelines'}

##########################################################
##########################################################
# Classifiers
classifiers = """
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

setup(
    # package information
    name='sc_pipelines',
    version=version,

    description='sc_pipelines : Dendrou group single-cell pipelines, master version',
    author='Charlotte Rich-Griffin',
    author_email='crg@well.ox.ac.uk',
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description='''sc_pipelines : Dendrou group single-cell pipelines, master version''',

    classifiers=[_f for _f in classifiers.split("\n") if _f],
    url="",
    # package contents
    packages=sc_pipelines_packages,
    package_dir=sc_pipelines_package_dirs,
    include_package_data=True,
    entry_points={
        "console_scripts": ["sc_pipelines = scpipelines.entry:main"]
    },
    # other options
    zip_safe=False,
    test_suite="tests",
)
