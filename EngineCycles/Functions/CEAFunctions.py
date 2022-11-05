from rocketcea.cea_obj import CEA_Obj
from rocketcea.cea_obj_w_units import CEA_Obj as CEA_Obj_w_units
import re
from typing import Optional
from functools import wraps


def get_values_from_cea_output(variable_name: str, column: int, full_output: str) -> float:
    match = re.search(fr'(?<={variable_name}) +([\s0-9.-]+)+', full_output)
    if match is None:
        raise ValueError(f'No match found in full output for [{variable_name}]')
    return float(match.group(0).split()[column])
    # return [float(group.replace(' ', '')) for group in match.groups()]


def cea_in_si_units(func):
    @wraps(func)
    def wrapper_func(**kwargs):
        kwargs['Pc'] /= 6894.757293  # Pascal to PSIA
        output = func(**kwargs)
        unit_conversion_tuple = (('c_star', 0.3048),  # feet to meter
                                 ('mm_cc', 1E-3),  # gram/mol to kilogram/mol
                                 ('mu_cc', 1E-4),  # milipoise to Pascal second
                                 ('cp_cc', 4184))  # From calorie/(gram Kelvin) to Joule/(kilogram Kelvin)
        for key, unit_factor in unit_conversion_tuple:
            try:
                output[key] *= unit_factor
            except KeyError:
                continue
        return output

    return wrapper_func


def cea_u_in_si_units(func):
    @wraps(func)
    def wrapper_func(**kwargs):
        kwargs['Pc'] *= 1E-5  # Pascal to Bar
        output = func(**kwargs)
        unit_conversion_tuple = (('mm_cc', 1E-3),  # gram/mol to kilogram/mol
                                 ('mu_cc', 1E-4),  # milipoise to Pascal second
                                 ('cp_cc', 1E+3))   # From kiloJoule/(kilogram Kelvin) to Joule/(kilogram Kelvin)
        for key, unit_factor in unit_conversion_tuple:
            try:
                output[key] *= unit_factor
            except KeyError:
                continue
        return output

    return wrapper_func


complete_regex_dict = {
    'c_star': ('CSTAR, M/SEC', 1),
    'C_F': (r'CF', 1),
    'T_C': ('T, K', 0),
    'mm_cc': ('M, [(]1/n[)]', 0),
    'y_cc': ('GAMMAs', 0),
    'mu_cc': ('VISC,MILLIPOISE', 0),
    'pr_cc': ('PRANDTL NUMBER', 0),
    'cp_cc': ('REACTIONS\n\n Cp, KJ/[(]KG[)][(]K[)]', 0)}


@cea_u_in_si_units
def get_cea_dict(fuelName: str, oxName: str, regex_dict: Optional[dict] = None, **kwargs):
    cea = CEA_Obj(fuelName=fuelName, oxName=oxName)
    full_output = cea.get_full_cea_output(**kwargs, short_output=1, pc_units='bar', output='siunits')
    if regex_dict is None:
        regex_dict = complete_regex_dict
    return {key: get_values_from_cea_output(variable_name=value[0],
                                            column=value[1],
                                            full_output=full_output)
            for key, value in regex_dict.items()}


def get_cea_chamber_dict(**kwargs):
    regex_dict = complete_regex_dict.copy()
    del regex_dict['C_F']
    del regex_dict['c_star']
    return get_cea_dict(regex_dict=regex_dict, **kwargs, eps=None)


def get_cea_performance_dict(**kwargs):
    regex_dict = {'c_star': ('CSTAR, FT/SEC', 1),
                  'C_F': (r'CF', 1)}
    return get_cea_dict(regex_dict=regex_dict, **kwargs)


if __name__ == '__main__':
    pass
    # kwargs1 = {'Pc': 55e5, 'MR': 5.6, 'eps': 22}
    # print(get_cea_dict(fuelName='LH2_NASA', oxName='LO2_NASA', **kwargs1))
    # get_cea_chamber_dict(fuelName='LH2_NASA', oxName='LO2_NASA', Pc=55e5, MR=5.6)
