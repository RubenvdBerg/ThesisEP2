from dataclasses import dataclass
from CoolProp.CoolProp import PropsSI


@dataclass
class A:
    @property
    def mass(self):
        return 1


@dataclass
class B(A):
    @property
    def mass(self):
        return 2


@dataclass
class C(B):
    @property
    def mass(self):
        return super(type(super()), self).mass


if __name__ == '__main__':
    h1 = PropsSI('H',
                 'T', 272.74,
                 'P', 15.5e6,
                 'n-Dodecane')
    h2 = PropsSI('H',
                 'T', 272.74 + 40,
                 'P', .25e6,
                 'n-Dodecane')
    Q = 1794693.7054650292
    m = 21.048244298120323
    dh = Q/m
    dh2 = h2 - h1
    print()
