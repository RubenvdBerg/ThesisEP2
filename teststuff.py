from dataclasses import dataclass

@dataclass
class Class1:
    attribute1: int
    attribute2: int = 2

@dataclass
class Class2(Class1):
    attribute3: int

if __name__ == '__main__':
    print(Class2(attribute1=1, attribute3=3))