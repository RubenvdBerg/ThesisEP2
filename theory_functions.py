# Liquid oxygen/liquid hydrogen (GH2 injection) 0.56-0.71
# Liquid oxygen/liquid hydrogen (LH2 injection) 0.76-1.02
# Liquid oxygen/RP-1 1.02-1.27

def average(iterable):
    return sum(iterable)/len(iterable)


if __name__ == '__main__':
    for values in [(0.56,0.71), (0.76, 1.02), (1.02, 1.27)]:
        print(average(values))
    options = {
    }
    print(options['hi'])