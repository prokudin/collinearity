"""
This module provides the interface to the configuration settings.

*conf* is a dictionary containing all of the shared information.
"""

conf = {}


def load_config(fname):
    """Load the settings stored in the file at *fname* into *conf*."""
    with open(fname) as f:
        L = f.readlines()
    D = {}
    for l in L:
        exec l in D
    conf.update(D['conf'])
