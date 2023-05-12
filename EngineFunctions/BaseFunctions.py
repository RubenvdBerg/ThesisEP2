from typing import Optional

import matplotlib.pyplot as plt
from numpy import log10


def only_one_none(a, b, c):
    a, b, c = a is None, b is None, c is None
    return a ^ b ^ c ^ all((a, b, c))


def multi_legend(axes: tuple[plt.Axes, ...], **kwargs):
    lines_list = []
    labels_list = []
    for ax in axes:
        lines, labels = ax.get_legend_handles_labels()
        for line, label in zip(lines, labels):
            lines_list.append(line)
            labels_list.append(label)
    axes[0].legend(lines_list, labels_list, **kwargs)


def copy_without(origin_dict, iterable_keys):
    copy_dict = origin_dict.copy()
    for key in iterable_keys:
        copy_dict.pop(key)
    return copy_dict


def format_si(value: float, unit: str, digits: int = 5, force_prefix: Optional[dict] = None):
    if force_prefix is None:
        force_prefix = {'g/s': 'k', 'Pa': 'M', 'W': 'M'}

    si_prefixes = ('Y', 'Z', 'E', 'P', 'T', 'G', 'M', 'k', '', 'm', '\u03BC', 'n', 'p', 'f', 'a', 'z', 'y')
    si_index = 8
    if value == 0:
        return 0
    n_before_comma = log10(abs(value)) // 1 + 1
    if force_prefix is not None and unit in force_prefix:
        prefix = force_prefix[unit]
        x = si_index - si_prefixes.index(prefix)
    else:
        x = round(n_before_comma / 3 - 1)
        if (1 > x > -2):
            x = 0
        si_index -= x
        prefix = si_prefixes[si_index]
    value *= 10 ** (-3 * x)
    decimals = int(digits - (n_before_comma - (3 * x)))
    format = min(digits, decimals)
    if digits == decimals:
        format = digits - 1
    val = f'{value:.{format}f}'
    if len(val) > digits + 1:
        val = f'{value:.{format - 1}f}'
    return rf'{val} {prefix}{unit}'


def format_fancy_name(name: str) -> str:
    return name.lower().replace(' ', '_')


def format_attr_name(attr_name: str) -> str:
    return attr_name.replace('_', ' ').replace('.', ' ').title()

def get_unit(attribute_name: str):
    unit_dict = {
        'energy_source_ratio': 'kg/s',
        'cc_prop_group_ratio': 'kg/s',
        '_mass': 'kg',
        '.mass': 'kg',
        'specific_impulse': 's',
        'pressure': 'Pa',
        'temp': 'K',
        'mass_flow': 'kg/s',
        'heat_capacity_ratio': '-',
        'specific_heat_capacity': 'J/kg/K',
        'molar_mass': 'kg/mol',
        'power': 'kg/s',
        'time': 's',
        'ratio': '-',
        'velocity': 'm/s',
        'delta_v': 'm/s',
    }
    for key, value in unit_dict.items():
        if key in attribute_name:
            return value


def get_symbol(attribute_name: str):
    attr_switcher = {
        'initial_mass': r'$m_0$',
        'overall_specific_impulse': r'$I_{sp}$',
        'power_mass': r'$m_{pow}$',
        'combustion_chamber_pressure': r'$p_{cc}$',
        'turbine_maximum_temperature': r'$(T_{tu})_{max}$',
        'burn_time': r'$t_b$',
        'change_in_velocity': r'$\Delta V$',
        'ideal_delta_v': r'$\Delta V$',
    }
    return attr_switcher[attribute_name]
