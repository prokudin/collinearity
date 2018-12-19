#!/usr/bin/env python
import sys
import os
import numpy as np
import time
import fnmatch
import cPickle
import zlib


def checkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def tex(x):
    return r'$\mathrm{' + x + '}$'


def save(data, name):
    compressed = zlib.compress(cPickle.dumps(data))
    with open(name, "wb") as f:
        f.writelines(compressed)


def load(name):
    with open(name, "rb") as compressed:
        data = cPickle.loads(zlib.decompress(compressed.read()))
    return data


def isnumeric(value):
    try:
        int(value)
        return True
    except:
        return False


def load_config(fname):
    with open(fname) as f:
        for l in f:
            exec l.replace('<<', '').replace('>>', '')
    return conf


def lprint(msg):
    sys.stdout.write('\r')
    sys.stdout.write(msg)
    sys.stdout.flush()
