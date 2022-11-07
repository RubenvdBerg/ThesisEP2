import os

from EngineCycles.ElectricPumpCycle.EPCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle.GGCycle import GasGeneratorCycle
from EngineCycles.BaseEngineCycle.EngineCycle import EngineCycle
from EngineCycles.OpenExpanderCycle.OECycle import OpenExpanderCycle
from KwakFix.KwakFixCycles import KwakFixElectricPumpCycle, KwakFixGasGeneratorCycle
from typing import Optional, Tuple, Callable
import arguments as args
import matplotlib.pyplot as plt


def get_engines(cycle_name: str, thrust: float = 100e3, burn_time: float = 300, p_cc_range: tuple = tuple(range(3, 11)),
                is_frozen: bool = False, verbose: bool = False, kwak: bool = False, _p_cc_factor: float = 1e6):
    extra_args_selection = {'ep': args.ep_arguments,
                            'gg': args.gg_arguments,
                            'oe': args.oe_arguments, }
    if kwak:
        cycle_selection = {'ep': KwakFixElectricPumpCycle,
                           'gg': KwakFixGasGeneratorCycle, }
        base_args = args.common_arguments_kwak
    else:
        cycle_selection = {'ep': ElectricPumpCycle,
                           'gg': GasGeneratorCycle,
                           'oe': OpenExpanderCycle, }
        base_args = args.base_arguments
    cycle = cycle_selection[cycle_name]
    extra_args = extra_args_selection[cycle_name]
    # noinspection PyArgumentList
    return [cycle(combustion_chamber_pressure=p_cc * _p_cc_factor,
                  thrust=thrust,
                  burn_time=burn_time,
                  is_frozen=is_frozen,
                  verbose=verbose,
                  **base_args, **extra_args) for p_cc in p_cc_range]


def ratio_plot(attribute_name: str, ylabel: str, title: str, pngname: str,
               engine_data: dict, p_cc_range: tuple, is_frozen: bool,
               line_val: Optional[float] = None,
               ylim: Optional[Tuple[float, float]] = None,
               sub_legend_kwargs: Optional[dict] = None,
               savefig: bool = False,
               dir_name: str = '', ):
    burn_times = list(engine_data)
    thrusts = list(engine_data[burn_times[0]])
    names = list(engine_data[burn_times[0]][thrusts[0]])

    fig, ax = plt.subplots()
    linestyles = ('-', '--', '-.', ':')
    colors = ('tab:red', 'tab:green', 'tab:orange', 'tab:blue', 'tab:purple')
    other_names = [name for name in names if name != 'EP']
    p = []

    for burn_time, color in zip(burn_times, colors):
        for thrust, linestyle in zip(thrusts, linestyles):
            label = f'$t_b${burn_time} $F_T${thrust}'
            for name in other_names:
                # Get the wanted attribute value for the EP-cycle and other cycle and find their ratio
                ep_values = [getattr(engine, attribute_name) for engine in engine_data[burn_time][thrust]['EP']]
                values = [getattr(engine, attribute_name) for engine in engine_data[burn_time][thrust][name]]
                ratios = [ep_val / val for ep_val, val in zip(ep_values, values)]
                # Plot attribute ratio and store legend info
                p_leg, = ax.plot(p_cc_range, ratios, label=label, color=color, linestyle=linestyle)
                p.append(p_leg)

    if ylim:
        ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('$p_{cc}$ [MPa]')
    mode = 'Frozen' if is_frozen else 'Shifting Equilibrium'
    ax.set_title(title + f'\n {mode}')
    if line_val is not None:
        ax.plot(p_cc_range, [line_val] * len(p_cc_range), color='black', linestyle='dotted')
    if sub_legend_kwargs is None:
        sub_legend_kwargs = {}
    vertical_header_legend(burn_times, thrusts, p, **sub_legend_kwargs)
    if savefig:
        plt.savefig(dir_name + '/' + pngname, dpi=1200)
    plt.show()


def vertical_header_legend(headervalues, subvalues, linelist,
                           dummy_pos=([3], [1]),
                           fontsize=8,
                           loc=0,
                           headerfunc=lambda x: fr'$t_b$:{x} s',
                           subfunc=lambda x: fr'$F_T$:{x} kN'):
    p5, = plt.plot(*dummy_pos, marker='None',
                   linestyle='None', label='dummy-tophead')
    categories = [subfunc(sub) for sub in subvalues]
    headerlabels = [headerfunc(header) for header in headervalues]
    labels = [item for headerlabel in headerlabels for item in (headerlabel, *categories)]
    l_sub, l_head = len(subvalues), len(headervalues)
    handles = [item for i in range(l_head) for item in [p5] + linelist[l_sub * i:l_sub * (i + 1)]]
    plt.legend(handles=handles,
               labels=labels,
               ncol=len(headervalues),
               loc=loc,
               fontsize=fontsize)


def plot_dv_ratio(default_ylim=True, **kwargs):
    if default_ylim:
        kwargs['ylim'] = (.72, .96)
    ratio_plot(attribute_name='ideal_delta_v',
               ylabel='$(I_{sp}*ln(MR^{-1}))_{ep}$ / $(I_{sp}*ln(MR^{-1}))_{gg}$ [-]',
               title='Comparison of ideal DeltaV of EP Cycle to GG Cycle',
               pngname='DeltaV_Comparison',
               sub_legend_kwargs={'fontsize': 7}, **kwargs)


def plot_mr_ratio(default_ylim=True, **kwargs):
    if default_ylim:
        kwargs['ylim'] = (1.2, 3.2)
    ratio_plot(attribute_name='mass_ratio',
               ylabel='$MR_{ep}/MR_{gg}$ [-]',
               title='Comparison of EP Cycle mass ratio to GG Cycle mass ratio',
               pngname='Mass_Ratio_Comparison',
               sub_legend_kwargs={'dummy_pos': ([3], [1.5]), 'loc': 2, 'fontsize': 6}, **kwargs)


def plot_m0_ratio(default_ylim=True, **kwargs):
    if default_ylim:
        kwargs['ylim'] = (.98, 1.010)
    ratio_plot(attribute_name='mass',
               ylabel=r'$m⁰_{ep}/m⁰_{gg}$ [-]',
               title='Initial mass comparison: EP Cycle to GG Cycle',
               pngname='Initial_Mass_Comparison',
               line_val=1, **kwargs)


def plot_isp_ratio(engine_data: dict, p_cc_range: tuple, is_frozen: bool, savefig=False, dir_name: str = '',
                   default_ylim: bool = True):
    colors = {'EP': 'tab:blue',
              'GG': 'tab:red',
              'OE': 'tab:green'}
    ratio_colors = colors.copy() | {'GG': 'tab:cyan'}
    markers = {'EP': 'o',
               'GG': 's',
               'OE': '^'}
    ratio_markers = {'GG': 'p',
                     'OE': 'H'}

    names = tuple(engine_data)
    other_names = [name for name in names if name != 'EP']
    isps = {}
    for name in names:
        isps[name] = [engine.overall_specific_impulse for engine in engine_data[name]]

    isp_ratios = {}
    for other_name in other_names:
        isp_ratios[other_name] = [ep_isp / isp for ep_isp, isp in zip(isps['EP'], isps[other_name])]

    fig, ax = plt.subplots()
    ax.set_xlabel('$p_{cc}$ [MPa]')
    ax.set_ylabel('$I_{sp}$ [s]')
    mode = 'Frozen' if is_frozen else 'Shifting Equilibrium'
    ax.set(title='Specific impulse comparison\n' + mode)
    ax2 = ax.twinx()
    ax2.set_ylabel('$I_{sp}$ Ratio (EP/GG)')
    if default_ylim:
        ax.set_ylim(340, 380)
        ax2.set_ylim(1.00, 1.08)

    lines = []
    for name in names:
        lines += ax.plot(p_cc_range, isps[name],
                         label=r'$I_{sp}$ ' + f'{name}-cycle',
                         color=colors[name],
                         marker=markers[name],
                         linestyle='dashed', )
        if name != 'EP':
            lines += ax2.plot(p_cc_range, isp_ratios[name],
                              label='$I_{sp}$ ratio ' + f'EP/{name}',
                              color=ratio_colors[name],
                              marker=ratio_markers[name],
                              markersize=8, )

    labels = [line.get_label() for line in lines]
    ax.legend(lines, labels, loc=2)
    if savefig:
        plt.savefig(dir_name + '/' + 'Specific_Impulse_Comparison', dpi=1200)
    plt.show()


def plot_all_ratio_plots(short: bool = False,
                         is_frozen: bool = False,
                         burn_times: tuple = (300, 600, 900, 1200),
                         thrusts: tuple = (25, 50, 75, 100),
                         names: tuple = ('EP', 'GG'),
                         p_cc_range: tuple = tuple(range(3, 11)),
                         dir_name: str = '',
                         savefig: bool = False,
                         default_ylim: bool = True,
                         **kwargs):
    if short:
        thrusts = [max(thrusts)]

    engine_data = {}
    for burn_time in burn_times:
        engine_data[burn_time] = {}
        for thrust in thrusts:
            engine_data[burn_time][thrust] = {}
            for name in names:
                engine_data[burn_time][thrust][name] = get_engines(cycle_name=name.lower(),
                                                                   thrust=thrust,
                                                                   burn_time=burn_time,
                                                                   is_frozen=is_frozen,
                                                                   p_cc_range=p_cc_range, **kwargs)
    other_kwargs = {'dir_name': dir_name,
                    'savefig': savefig,
                    'default_ylim': default_ylim,
                    'is_frozen': is_frozen,
                    'p_cc_range': p_cc_range}

    plot_isp_ratio(engine_data=engine_data[max(burn_times)][max(thrusts)], **other_kwargs)
    plot_dv_ratio(engine_data=engine_data, **other_kwargs)
    plot_mr_ratio(engine_data=engine_data, **other_kwargs)
    plot_m0_ratio(engine_data=engine_data, **other_kwargs)


if __name__ == '__main__':
    kwargs = {'p_cc_range': tuple(range(3,11)),
              'is_frozen': False,
              'kwak': False}
    engine_data = {}
    for name in ('EP', 'GG'):
        engine_data[name] = get_engines(cycle_name=name.lower(), **kwargs)
    plot_isp_ratio(engine_data=engine_data, p_cc_range=kwargs['p_cc_range'], is_frozen=kwargs['is_frozen'] )
    # plot_all_ratio_plots(short=True, kwak=True)
