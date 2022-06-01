from dataclasses import dataclass
# def array_create(n):
#     return [n + 1, n + 2], [n + 3]
#
#
# long_array = [array_create(n) for n in range(5)]
# print([f[0] for x, f in long_array])

@dataclass
class One:
    one: float
    def __post_init__(self):
        self.two = 2

@dataclass
class Two(One):
    three: float
    def __post_init__(self):
        self.four = 4
        super().__post_init__()

if __name__ == '__main__':
    two = Two(one=1, three=3)
    print(two.two)