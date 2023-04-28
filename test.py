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
    # h1 = PropsSI('H',
    #              'T', 272.74,
    #              'P', 15.5e6,
    #              'n-Dodecane')
    # h2 = PropsSI('H',
    #              'T', 272.74 + 40,
    #              'P', .25e6,
    #              'n-Dodecane')
    # Q = 1794693.7054650292
    # m = 21.048244298120323
    # dh = Q/m
    # dh2 = h2 - h1
    # print()
    a = [1,2,3]
    b = a + [
        4,
        5,
        6,
    ]
    print(b)
    temp, pressure = 922, 4.7e6
    h2 = PropsSI('CPMOLAR',
                'T', temp,
                'P', pressure,
                'Hydrogen')
    h2o = PropsSI('CPMOLAR',
                'T', temp,
                'P', pressure,
                'Water')
    cp = .12*h2o+.88*h2
    mix = PropsSI('M',
                'T', temp,
                'P', pressure,
                'Water[.12]&Hydrogen[.88]')
    print(cp, mix, cp/mix)
