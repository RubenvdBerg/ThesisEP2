from dataclasses import dataclass


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
    c = C()
    print(c.mass)
