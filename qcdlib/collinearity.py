# -*- coding: utf-8 -*-
from weakref import WeakValueDictionary
import numpy as np

DEFAULTS = {}
DEFAULTS["hadron"] = "pi+"
DEFAULTS["M"] = 0.938
DEFAULTS["MJ"] = 0.3
DEFAULTS["MX"] = 1.3
DEFAULTS["Ma"] = 1.5
DEFAULTS["Mb"] = 0.3
DEFAULTS["MfT2"] = 0.5
DEFAULTS["MiT2"] = 0.5
DEFAULTS["deltaM"] = 0.3
DEFAULTS["delta_kT"] = 0.3
DEFAULTS["kf"] = 0.3
DEFAULTS["ki"] = 0.3
DEFAULTS["kT"] = 0.3

_collinearity_cache = WeakValueDictionary()


class _Collinearity(object):
    def __init__(self, **kwargs):
        diff = sorted(kwargs.keys() - DEFAULTS.keys())

        if diff:
            raise ValueError("Invalid Keyword: '{0}'".format(diff[0]))

        for k in DEFAULTS:
            if k in kwargs:
                setattr(self, k, kwargs[k])
            else:
                setattr(self, k, DEFAULTS[k])

        self.kwargs = sorted(kwargs.items())

        self.kf2 = self.kf ** 2
        self.ki2 = self.ki ** 2
        self.kT2 = self.kT ** 2

        self.M2 = self.M ** 2

        self.MfT = self.MfT2 ** 0.5
        # self.MfT = np.sqrt(self.kT ** 2 + self.MJ ** 2 + self.deltaM ** 2)
        # XXX: ^^^ which to use?

        self.MiT = self.MiT2 ** 0.5

        self.delta_kT2 = self.delta_kT ** 2

        if self.hadron in ('pi+', 'pi-', 'pi0'):
            self.Mh = 0.135
        elif self.hadron in ('k+', 'k-', 'k0'):
            self.Mh = 0.493
        else:
            raise ValueError("hadron")

        self.Mh2 = self.Mh ** 2

    def __repr__(self):
        return "Collinearity(%s)" % ", ".join("%s=%s" % (str(k), repr(v))
                                              for k, v in self.kwargs)

    def get_W2(self, x, Q2):
        W2 = Q2 * (1. - x) / x + self.M2
        return W2

    def get_MBT(self, PhT):
        MBT = np.sqrt(self.Mh2 + PhT ** 2)
        return MBT

    def get_MiT(self, x, Q2):
        xn = self.get_xn(x, Q2)
        MiT = np.sqrt((xn * self.kT**2 + xn *
                       (self.Ma + self.Mb / np.sqrt(xn)) ** 2 -
                       (1 - xn) * xn * self.M2 + self.deltaM ** 2) / (1 - xn))
        return MiT

    def get_xn(self, x, Q2):
        return 2 * x / (1 + np.sqrt(1 + 4 * x**2 * self.M2 / Q2))

    def get_yh(self, x, z, Q2, PhT, *, sign=-1):
        xn = self.get_xn(x, Q2)
        expy = (Q2**0.5 * z * (Q2 - xn**2 * self.M2) /
                (2 * self.M2 * xn**2 * (self.Mh2 + PhT**2)**0.5) +
                sign * Q2**0.5 / (xn * self.M) *
                (z**2 * (Q2 - xn**2 * self.M2)**2 /
                (4 * self.M2 * xn**2 * (self.Mh2 + PhT**2)) - 1)**0.5)
        return np.log(expy)

    def get_zn(self, x, z, Q2, PhT):
        xn = self.get_xn(x, Q2)
        MBT = self.get_MBT(PhT)
        zn = (xn * z / (2 * x) * (1. +
              np.sqrt(1. - 4. * self.M2 * MBT**2 * x**2 / (Q2**2 * z**2))))
        return zn

    def get_yp(self, x, Q2):
        xn = self.get_xn(x, Q2)
        return np.log(np.sqrt(Q2) / (xn * self.M))

    def get_yi(self, Q2):
        return 0.5 * np.log(Q2 / self.MiT**2)

    def get_yf(self, Q2):
        return -0.5 * np.log(Q2 / self.MfT**2)

    def get_MhT(self, PhT):
        return np.sqrt(self.Mh2 + PhT**2)

    def get_R(self, x, z, Q2, PhT):
        MiT = self.get_MiT(x, Q2)
        yi = self.get_yi(Q2)
        yf = self.get_yf(Q2)
        MhT = self.get_MhT(PhT)
        yh = self.get_yh(x, z, Q2, PhT)
        zn = self.get_zn(x, z, Q2, PhT)
        znhat = zn / z
        Ph_kf = (0.5 * MhT * self.MfT * (np.exp(yf - yh) + np.exp(yh - yf)) -
                 znhat / zn * PhT**2 - PhT * self.kT)
        Ph_ki = (0.5 * MhT * MiT * (np.exp(yi - yh) - np.exp(yh - yi)) -
                 PhT * self.kT)
        return np.abs(Ph_kf / Ph_ki)

    # R0, Eq. (4.14)
    def get_R0(self, Q2):
        """Collinearity ratio defined in the paper Eq. (4.15)"""
        return max(self.ki2 / Q2, self.kf2 / Q2, self.kT2 / Q2)

    # We call R1 in the new paper what was R in the previous ...
    def get_R1(self, x, z, Q2, PhT):
        """Collinearity ratio defined in the paper Eq. (4.15)"""
        return self.get_R(x, z, Q2, PhT)

    # R2, Eq. (4.17)
    def get_R2(self, x, z, Q2, PhT):
        zn = self.get_zn(x, z, Q2, PhT)
        znhat = zn / z
        qT = -PhT / zn
        return np.abs(-(1. - znhat) - znhat * qT**2 / Q2 - (1. - znhat) *
                      self.kf2 / (znhat * Q2) - self.delta_kT2 /
                      (znhat * Q2) + 2. * qT * self.delta_kT / Q2)


def Collinearity(**kwargs):
    kw_hash = tuple(sorted(kwargs.items()))

    if kw_hash not in _collinearity_cache:
        c = _Collinearity(**kwargs)
        _collinearity_cache[kw_hash] = c
    else:
        c = _collinearity_cache[kw_hash]

    return c


if __name__ == '__main__':
    kin = Collinearity(hadron="pi+")

    x = 0.01
    z = 0.5
    Q2 = 4.0
    PhT = 0.1

    R0 = kin.get_R0(Q2)
    R1 = kin.get_R1(x, z, Q2, PhT)
    R2 = kin.get_R2(x, z, Q2, PhT)

    print("R0:", R0)
    print("R1:", R1)
    print("R2:", R2)
