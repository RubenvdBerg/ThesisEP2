def shomate_equation(temperature, a, b, c, d, e):
    t = temperature / 1000
    return a + b * t + c * t ** 2 + d * t ** 3 + e / t ** 2


def kdb_equation(temp, a, b, c, d, e):
    return a + b * temp + c * temp ** 2 + d * temp ** 3 + e * temp ** 4


def get_data(name, temp):
    coeff_dict = {}
    for key in data[name]:
        if key[0] <= temp <= key[1]:
            coeff_dict = data[name][key]
    if coeff_dict == {}:
        raise ValueError(f'{temp} is out of range for {name}, no coefficients could be found')
    elif all(item is None for item in coeff_dict.items()):
        raise ValueError(f'No coefficient data is available for {name}')
    return coeff_dict


def get_heat_capacity(name, temp):
    coeff_dict = get_data(name=name, temp=temp)
    return kdb_equation(temp, **coeff_dict)


# Input in [K], output in [J/(mol*K)]
data = {'LO2': {(54.75, 143.15): {'a': -6.141212E+01,
                                  'b': 4.329776E+00,
                                  'c': -5.289717E-02,
                                  'd': 2.109594E-04,
                                  'e': 0.}},
        'GO2': {(0., 10000000.): {'a': 2.970450E+01,
                                  'b': -9.895231E-03,
                                  'c': 3.989792E-05,
                                  'd': -3.394227E-08,
                                  'e': 9.184016E-12}},
        'LH2': {(13.75, 38.150): {'a': 3.196527E+01,
                                  'b': -2.781568E+00,
                                  'c': 1.026507E-01,
                                  'd': -2.053536E-05,
                                  'e': 0.}},
        'GH2': {(0., 10000000.): {'a': 2.700357E+01,
                                  'b': 1.193388E-02,
                                  'c': -2.407279E-05,
                                  'd': 2.146124E-08,
                                  'e': -6.147960E-12}},
        'LRP1': {(260.00, 690.0): {'a': -220.5821,
                                   'b': 4.218693E0,
                                   'c': -1.219061E-2,
                                   'd': 1.723811E-5,
                                   'e': -8.296845E-9}},
        'GRP1': {(0., 10000000.): {'a': None,
                                   'b': None,
                                   'c': None,
                                   'd': None,
                                   'e': None}},
        }
