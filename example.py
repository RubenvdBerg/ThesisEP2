class A:
    variable_a: float


class B1(A):
    variable_b1: int

    def method_b1(self):
        pass


class B2(A):
    variable_b2: int

    def method_b2(self):
        pass


class C:
    variable_c: int

    def method_c(self):
       pass


class C1(C, B1):
    pass


class C2(C, B2):
    pass

