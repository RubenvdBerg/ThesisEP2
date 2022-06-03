from dataclasses import dataclass, field
from typing import Optional


def get_attribute2(attribute1: float):
    return attribute1 * 100


def get_attribute1(attribute2: float):
    return attribute2 / 100


@dataclass
class Class1:
    attribute1: Optional[float] = None
    attribute2: Optional[float] = None

    def __post_init__(self):
        if not (self.attribute1 is None) ^ (self.attribute2 is None):
            raise ValueError(
                'Neither or both attribute1 and attribute2 are given. Provide one and only one')
        elif self.attribute1 is None:
            self.attribute1 = get_attribute1(self.attribute2)
        elif self.attribute2 is None:
            self.attribute2 = get_attribute2(self.attribute1)


@dataclass
class Class2:
    attribute1_input: Optional[float] = None
    attribute2_input: Optional[float] = None

    def __post_init__(self):
        if not (self.attribute1_input is None) ^ (self.attribute2_input is None):
            raise ValueError(
                'Neither or both attribute1 and attribute2 are given. Provide one and only one')

    @property
    def attribute1(self):
        if self.attribute1_input is None:
            return get_attribute1(self.attribute2)
        else:
            return self.attribute1_input

    @property
    def attribute2(self):
        if self.attribute2_input is None:
            return get_attribute2(self.attribute1)
        else:
            return self.attribute2_input

    @attribute1.setter
    def attribute1(self, attribute1: float):
        self.attribute2_input = get_attribute2(attribute1)
        self.attribute1_input = attribute1

    @attribute2.setter
    def attribute2(self, attribute2: float):
        self.attribute1_input = get_attribute1(attribute2)
        self.attribute2_input = attribute2


if __name__ == '__main__':
    ec1 = Class1(attribute1=1.)
    ec2 = Class2(attribute1_input=1.)

    for ec in (ec1, ec2):
        ec.attribute2 = 200
        print(f'{ec.attribute1} should be 2.0')
