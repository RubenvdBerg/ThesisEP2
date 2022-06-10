from dataclasses import dataclass


@dataclass
class Class1:
    attr1: float = 1
    attr2: float = 2
    attr3: float = 3

    def __post_init__(self):
        self.attr10 = 10

    @property
    def prop1(self):
        return 1

@dataclass
class Class2(Class1):
    attr4: bool = True

    def __post_init__(self):
        super().__post_init__()
        self.attr11 = 11

    @property
    def prop2(self):
        return 2

@dataclass
class Class4(Class1):
    attr5: float = 0

    @property
    def prop1(self):
        return 3


@dataclass
class Class3(Class4, Class2):
    attr6: float = 6
    pass


def make_error_message(attribute_name, object):
    attribute_value = getattr(object, attribute_name)
    return f'{attribute_name} is already given with value {attribute_value}'


if __name__ == '__main__':
    c = Class1
    print(make_error_message('attr1', c))
    print(Class3(attr4=False).attr10)
