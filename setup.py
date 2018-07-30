"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
import re

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Read the file version (see http://stackoverflow.com/a/7071358)
PACKAGE_NAME = "pyp_beagle"
PACKAGE_DIR = "PyP-BEAGLE"
VERSION_FILE = PACKAGE_DIR + "/_version.py"
verstrline = open(VERSION_FILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSION_FILE,))

setup(
    name=PACKAGE_NAME,
    version=verstr,
    description='Package for post-processing of results obtained with the Beagle SED fitting tool',
    long_description=long_description,
    url='https://github.com/jacopo-chevallard/PyP-BEAGLE',
    author='Jacopo Chevallard',
    author_email='jacopo.chevallard@gmail.com',
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7',
    ],

    # What does your project relate to?
    keywords='astronomy galaxies statistics visualization',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=[PACKAGE_NAME],

    package_dir={PACKAGE_NAME: PACKAGE_DIR},

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["postprocess_beagle_results.py"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['matplotlib', 'scipy', 'numpy', 'atpy', 'getdist', 'pathos', 'astropy', 'bokeh', 'natsort'],

    include_package_data=True, 

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            PACKAGE_NAME+'='+PACKAGE_NAME+'.command_line:main',
        ],
    },
)
