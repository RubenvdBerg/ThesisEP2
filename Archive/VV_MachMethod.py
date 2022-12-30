from math import sqrt
def get_initial_estimate_mach(local_area_ratio, is_subsonic=False, heat_capacity_ratio=1.14):
    a_at = local_area_ratio
    y = heat_capacity_ratio
    p = 2 / (y + 1)
    q = 1 - p
    if is_subsonic:
        r, a, s = a_at ** 2, p ** (1 / q), 1
    else:
        r, a, s = a_at ** (2 * q / p), q ** (1 / p), -1
    r2 = (r - 1) / (2 * a)
    initial_guess = 1 / ((1 + r2) + sqrt(r2 * (r2 + 2)))

    def get_f_and_derivs(x):
        f = (p + q * x) ** (1 / q) - r * x
        df = (p + q * x) ** (1 / q - 1) - r
        ddf = p * ((p + q * x) ** (1 / q - 2))
        return f, df, ddf

    def newton_raphson_plus(x):
        while True:
            f, df, ddf = get_f_and_derivs(x)
            xnew = x - 2 * f / (df - sqrt(df**2 - 2 * f * ddf))
            if abs(xnew - x) / xnew < .0001:
                break
            x = xnew
        return xnew
    final = newton_raphson_plus(initial_guess)
    return sqrt(final)**s

if __name__ == '__main__':
    print(get_initial_estimate_mach(local_area_ratio=2.085, is_subsonic=True))