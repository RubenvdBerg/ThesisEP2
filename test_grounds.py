# from dataclasses import dataclass, field
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
a =1
b = -a
print(b)