def array_create(n):
    return [n + 1, n + 2], [n + 3]


long_array = [array_create(n) for n in range(5)]
print([f[0] for x, f in long_array])
