import numpy as np
from rocketcea.cea_obj import CEA_Obj
from rocketcea.cea_obj_w_units import CEA_Obj as CEA_Obj_w_units
import re
from typing import Optional
from functools import wraps
from numpy import logspace, interp
from EngineFunctions.EmpiricalRelations import get_gas_generator_mmr_rp1


def get_match_from_cea_output(variable_name: str, full_output: str) -> list:
    match = re.findall(fr'(?<={variable_name}) +([\s0-9.-]+)+', full_output)
    if match is None:
        raise ValueError(f'No match found in full output for [{variable_name}]')
    return match


def get_values_from_cea_output(variable_name: str, column: int, full_output: str) -> float:
    match = get_match_from_cea_output(variable_name=variable_name, full_output=full_output)
    values = re.findall(r'[0-9.]+', match[0])
    exponents = re.findall(r'[ -][0-9]+ ', match[0])
    value = float(values[column])
    if exponents:
        exponent = float(exponents[column])
        return value * 10 ** exponent
    else:
        return value


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
    return {'full_output': full_output} | {key: get_values_from_cea_output(variable_name=value[0],
                                                                           column=value[1],
                                                                           full_output=full_output)
                                           for key, value in regex_dict.items()}


def get_cea_dict_gg(**kwargs):
    return get_cea_dict(regex_dict={'y_cc': ('GAMMAs', 0),
                                    'cp_cc': ('REACTIONS\n\n Cp, KJ/[(]KG[)][(]K[)]', 0),
                                    'mm_cc': ('M, [(]1/n[)]', 0),
                                    # 'mm_cc2': ('MW, MOL WT', 0),
                                    'T_C': ('T, K', 0),
                                    'rho_cc': ('RHO, KG/CU M', 0)},

                        eps=None,
                        **kwargs)


def get_gas_generator_mmr(temperature_limit: float, fuelName: str, oxName: str, Pc: float):
    if 'LH2' in fuelName:
        range_tuple = (.01, 6.)
    elif 'RP' in fuelName:
        range_tuple = (.01, 3.)
    elif 'CH4' in fuelName:
        range_tuple = (.01, 4.)
    cea_obj = CEA_Obj(fuelName=fuelName, oxName=oxName)

    def get_t_comb(MR: float, Pc: float):
        Pc /= 6894.76  # Pa to PSIA
        t_comb_rankine = cea_obj.get_Tcomb(Pc=Pc, MR=MR)
        return t_comb_rankine / 1.8  # Rankine to Kelvin

    mixture_ratios = np.linspace(*range_tuple, num=100)
    combustion_temperatures = [get_t_comb(MR=mr, Pc=Pc) for mr in mixture_ratios]
    cea_mmr = np.interp(temperature_limit, combustion_temperatures, mixture_ratios)
    return cea_mmr



def get_cea_chamber_dict(**kwargs):
    regex_dict = complete_regex_dict.copy()
    del regex_dict['C_F']
    del regex_dict['c_star']
    return get_cea_dict(regex_dict=regex_dict, **kwargs, eps=None)


