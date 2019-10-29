# -*- coding: utf-8 -*-
from collections import namedtuple
import numpy as np

defaults = {
    "M": 0.938,
    "MJ": 0.3,
    "MX": 1.3,
    "Ma": 1.5,
    "Mb": 0.3,
    "MfT2": 0.5,
    "MiT2": 0.5,
    "deltaM": 0.3,
    "delta_kT": 0.3,
    "hadron": "pi+",
    "kT": 0.3,
    "kf": 0.3,
    "ki": 0.3
}

derived = {
    "M2": lambda x: x["M"] ** 2,
    "MfT": lambda x: x["MfT2"] ** 0.5,
    "Mh2": lambda x: (0.135 ** 2 if x["hadron"] in ("pi+", "pi-", "pi0") else
                      0.493 ** 2 if x["hadron"] in ("k+", "k-", "k0") else
                      None),
    "MiT": lambda x: x["MiT2"] ** 0.5,
    "delta_kT2": lambda x: x["delta_kT"] ** 2,
    "kT2": lambda x: x["kT"] ** 2,
    "kf2": lambda x: x["kf"] ** 2,
    "ki2": lambda x: x["ki"] ** 2
}

# MfT = np.sqrt(kT ** 2 + MJ ** 2 + deltaM ** 2)  # XXX: use which?

_Collinearity = namedtuple("_Collinearity",
                           sorted(defaults.keys() | derived.keys()))


class Collinearity(_Collinearity):
    __slots__ = ()

    def __new__(cls, **kwargs):
        diff = kwargs.keys() - defaults.keys()

        if diff:
            raise ValueError("Invalid Keyword: '{0}'".format(sorted(diff)[0]))

        d = defaults.copy()
        d.update(kwargs)

        for key, func in derived.items():
            d[key] = func(d)

        return super().__new__(Collinearity, **d)

    def get_MBT(self, PhT):
        MBT = np.sqrt(self.Mh2 + PhT ** 2)
        return MBT

    def get_MhT(self, PhT):
        return np.sqrt(self.Mh2 + PhT**2)

    def get_MiT(self, x, Q2):
        xn = self.get_xn(x, Q2)
        MiT = np.sqrt((xn * self.kT**2 + xn *
                       (self.Ma + self.Mb / np.sqrt(xn)) ** 2 -
                       (1 - xn) * xn * self.M2 + self.deltaM ** 2) / (1 - xn))
        return MiT

    def get_W2(self, x, Q2):
        W2 = Q2 * (1. - x) / x + self.M2
        return W2

    def get_xn(self, x, Q2):
        return 2 * x / (1 + np.sqrt(1 + 4 * x**2 * self.M2 / Q2))

    def get_yf(self, Q2):
        return -0.5 * np.log(Q2 / self.MfT**2)

    def get_yh(self, x, z, Q2, PhT, *, sign=-1):
        xn = self.get_xn(x, Q2)
        expy = (Q2**0.5 * z * (Q2 - xn**2 * self.M2) /
                (2 * self.M2 * xn**2 * (self.Mh2 + PhT**2)**0.5) +
                sign * Q2**0.5 / (xn * self.M) *
                (z**2 * (Q2 - xn**2 * self.M2)**2 /
                (4 * self.M2 * xn**2 * (self.Mh2 + PhT**2)) - 1)**0.5)
        return np.log(expy)

    def get_yi(self, Q2):
        return 0.5 * np.log(Q2 / self.MiT**2)

    def get_yp(self, x, Q2):
        xn = self.get_xn(x, Q2)
        return np.log(np.sqrt(Q2) / (xn * self.M))

    def get_zn(self, x, z, Q2, PhT):
        xn = self.get_xn(x, Q2)
        MBT = self.get_MBT(PhT)
        zn = (xn * z / (2 * x) * (1. +
              np.sqrt(1. - 4. * self.M2 * MBT**2 * x**2 / (Q2**2 * z**2))))
        return zn

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
