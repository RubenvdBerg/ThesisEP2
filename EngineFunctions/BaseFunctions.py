import warnings
from typing import Optional

import matplotlib.pyplot as plt
from numpy import log10


def only_one_none(a, b, c):
    a, b, c = a is None, b is None, c is None
    return a ^ b ^ c ^ all((a, b, c))


def multi_legend(axes: tuple[plt.Axes,...], **kwargs):
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

    si_prefixes = ('Y', 'Z', 'E', 'P', 'T', 'G', 'M', 'k', '', 'd', 'c', 'm', '\u03BC', 'n', 'p', 'f', 'a', 'z', 'y')
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