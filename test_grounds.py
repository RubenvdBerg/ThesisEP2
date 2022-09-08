from dataclasses import dataclass, field


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

    def __post_init__(self):
        setattr(self, 'attr2', 2)

    @property
    def whatever(self):
        return self.attr1 * 2

print(f'{1:{""}}')
d1 = Class1(attr1=1)
d1.attr1 = 2
print(d1.whatever)
print(d1.attr2)
print((1, 2) + (1,))
