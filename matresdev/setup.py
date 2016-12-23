import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="matresdev",
    version="0.0.4",
    author="Rostislav Chudoba",
    author_email="rostislav.chudoba@rwth-aachen.de",
    description=("Package supporting material research and development "
                 "with the focus on material characterization."),
    license="BSD",
    keywords="simple material database",
    url="http://packages.python.org/matresdev",
    packages=['matresdev'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
)
