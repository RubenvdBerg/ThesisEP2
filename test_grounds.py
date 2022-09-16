from dataclasses import dataclass, field
from typing import ClassVar

#
# @dataclass
# class Class1:
#     attr1: float = 1
#
# @dataclass
# class Class2:
#     attr2: Class1 = field(init=False, default=Class1())
#
# def default_factory():
#     return Class1()
#
# @dataclass
# class Class2_better:
#     attr2: Class1 = field(init=False, default_factory=default_factory)
#
# if __name__ == '__main__':
#     # c2 = Class2()
#     # c2.attr2.attr1 = 3
#     # c3 = Class2()
#     # print(c3.attr2.attr1)
# 
#     c2 = Class2_better()
#     c2.attr2.attr1 = 3
#     c3 = Class2_better()
#     print(c3.attr2.attr1)
@dataclass
class Class1:
    attr1: float
    clasvar: ClassVar[bool] = False

    def __post_init__(self):
        setattr(self, 'attr2', 2)
        Class1.clasvar = True

    @property
    def whatever(self):
        return self.attr1 * 2

@dataclass
class Class2:

    def __post_init__(self):
        self.class3._attr1 = 1

    @property
    def class3(self):
        return Class3()

@dataclass
class Class3:
    _attr1: float = 0

    @property
    def attr1(self):
        return self._attr1

c2 = Class2()
print(c2.class3.attr1)
print(f'{1:{""}}')
print(Class1.clasvar)
d1 = Class1(attr1=1)
print(Class1.clasvar)
d1.attr1 = 2
print(d1.whatever)
print(d1.attr2)
print((1, 2) + (1,))
