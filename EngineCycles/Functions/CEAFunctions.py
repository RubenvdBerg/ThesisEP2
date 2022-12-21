from rocketcea.cea_obj import CEA_Obj
from rocketcea.cea_obj_w_units import CEA_Obj as CEA_Obj_w_units
import re
from typing import Optional
from functools import wraps
from numpy import logspace, interp


def get_values_from_cea_output(variable_name: str, column: int, full_output: str) -> float:
    match = get_match_from_cea_output(variable_name=variable_name, full_output=full_output)
    return float(match[0].split()[column])


def get_match_from_cea_output(variable_name: str, full_output: str) -> list:
    match = re.findall(fr'(?<={variable_name}) +([\s0-9.-]+)+', full_output)
    if match is None:
        raise ValueError(f'No match found in full output for [{variable_name}]')
    return match


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
                                 ('cp_cc', 1E+3))  # From kiloJoule/(kilogram Kelvin) to Joule/(kilogram Kelvin)
        for key, unit_factor in unit_conversion_tuple:
            try:
                if type(output[key]) == float:
                    output[key] *= unit_factor
                else:
                    print()
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


def get_cea_dict_gg(**kwargs):
    return get_cea_dict(regex_dict={'y_cc': ('GAMMAs', 0),
                                    'cp_cc': ('REACTIONS\n\n Cp, KJ/[(]KG[)][(]K[)]', 0),
                                    'mm_cc': ('M, [(]1/n[)]', 0),
                                    'T_C': ('T, K', 0), },

                        eps=None,
                        **kwargs)


def get_gas_generator_mmr(fuelName: str, oxName: str, temp_limit: float, **kwargs):
    if 'LH2' in fuelName:
        log_range = (-.8, .4)
    elif 'RP1' in fuelName:
        log_range = (-1.60, .1)
    elif 'CH4' in fuelName:
        log_range = (-1.4, .15)
    mmr_list = list(float(f'{x:.2f}') for x in logspace(*log_range, 25))

    i = 0
    while True:
        # Increase the start of MMR range until CEA doesn't crash, then return temp range and find MMR from interpolation
        cea = CEA_Obj(fuelName=fuelName, oxName=oxName)
        kwargs['Pc'] *= 1e-5
        full_output = cea.get_full_cea_output(MR=mmr_list, **kwargs, short_output=1, pc_units='bar', output='siunits')
        match = get_match_from_cea_output(variable_name='T, K', full_output=full_output)
        temps = [float(val.split()[0]) for val in match]
        if temps:
            break

        i += 1
        mmr_list.pop(0)
        if i > 10:
            raise ValueError(
                'Could not solve gas generator mass mixture ratio (MMR) with CEA for given turbine temperature limit. Give the gas generator MMR manually or change the turbine temp limit.')
    # Ugly bug fix, sometimes full output contains a duplicate at the start
    if len(temps) > len(mmr_list):
        temps = temps[-len(mmr_list):]
    return interp(temp_limit, temps, mmr_list)


def get_cea_chamber_dict(**kwargs):
    regex_dict = complete_regex_dict.copy()
    del regex_dict['C_F']
    del regex_dict['c_star']
    return get_cea_dict(regex_dict=regex_dict, **kwargs, eps=None)


if __name__ == '__main__':
    kwargs1 = {'Pc': 34.5, 'eps': None, 'fuelName': 'CH4', 'oxName': 'LO2_NASA'}
    kwargs2 = kwargs1 | {'MR': (.32, .33, .34), 'PcOvPe': 100, }
    # print(get_cea_dict(fuelName='LH2_NASA', oxName='LO2_NASA', regex_dict={'T_C': ('T, K', 0)}, frozen=1, frozenAtThroat=1, **kwargs1))
    get_gas_generator_mmr(**kwargs1, frozen=1, frozenAtThroat=1, temp_limit=900)

    # get_cea_chamber_dict(fuelName='LH2_NASA', oxName='LO2_NASA', Pc=55e5, MR=5.6)
