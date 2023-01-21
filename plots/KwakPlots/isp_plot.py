import os

from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from EngineCycles.CoolantBleedCycle import CoolantBleedCycle
from KwakFix.KwakFixCycles import KwakFixElectricPumpCycle, KwakFixGasGeneratorCycle
from typing import Optional, Tuple
from EngineArguments import arguments as args
from EngineArguments.default_arguments import get_default_kwargs
import matplotlib.pyplot as plt
import matplotlib as mpl
from time import strftime
import json
import scienceplots


# plt.style.use('science')
# mpl.use('Qt5Agg')

def get_engines(cycle_name: str, thrust: float = 100e3, burn_time: float = 300, p_cc_range: tuple = tuple(range(3, 11)),
                is_frozen: bool = False, verbose: bool = False, kwak: bool = False, _p_cc_factor: float = 1e6,
                **kwargs):
    if kwak:
        cycle_selection = {'ep': KwakFixElectricPumpCycle,
                           'gg': KwakFixGasGeneratorCycle, }

    else:
        cycle_selection = {'ep': ElectricPumpCycle,
                           'gg': GasGeneratorCycle,
                           'oe': OpenExpanderCycle,
                           'cb': CoolantBleedCycle, }

    cycle = cycle_selection[cycle_name]

    default_args = get_default_kwargs(cycle)
    total_args = default_args | kwargs | {'thrust': thrust,
                                          'burn_time': burn_time,
                                          'is_frozen': is_frozen,
                                          'verbose': verbose, }

    def yield_engines():
        for p_cc in p_cc_range:
            # noinspection PyArgumentList
            yield cycle(combustion_chamber_pressure=p_cc * _p_cc_factor, **total_args)

    return tuple(yield_engines())


def get_data_dict(attributes: tuple, thrusts: tuple, burn_times: tuple, names: tuple,
                  p_cc_range: tuple = tuple(range(3, 11)), is_frozen: bool = False, kwak: bool = False, **kwargs):
    engine_data = {}
    for burn_time in burn_times:
        engine_data[burn_time] = {}
        for thrust in thrusts:
            engine_data[burn_time][thrust] = {}
            for name in names:
                engine_data[burn_time][thrust][name] = get_engines(cycle_name=name.lower(),
                                                                   thrust=thrust,
                                                                   burn_time=burn_time,
                                                                   p_cc_range=p_cc_range,
                                                                   is_frozen=is_frozen,
                                                                   kwak=kwak,
                                                                   **kwargs)
    data_dict = {}
    for attribute in attributes:
        data_dict[attribute] = {}
        for name in names:
            data_dict[attribute][name] = {}
            for burn_time in burn_times:
                data_dict[attribute][name][burn_time] = {}
                for thrust in thrusts:
                    data_dict[attribute][name][burn_time][thrust] = [getattr(engine, attribute) for engine in
                                                                     engine_data[burn_time][thrust][name]]

    data_dict['info'] = kwargs | {'burn_times': burn_times, 'thrusts': thrusts, 'names': names, 'is_frozen': is_frozen,
                                  'p_cc_range': p_cc_range, 'kwak': kwak}

    return data_dict


def ratio_plot(attribute_name: str, ylabel: str, variable_title: str, pngname: str, p_cc_range: tuple, is_frozen: bool,
               kwak: bool,
               data_dict: Optional[dict] = None,
               line_val: Optional[float] = None,
               ylim: Optional[Tuple[float, float]] = None,
               sub_legend_kwargs: Optional[dict] = None,
               savefig: bool = False,
               dir_name: str = '',
               unit: str = '-',
               ):
    linestyles = ('-', '--', '-.', ':')
    colors = ('tab:red', 'tab:green', 'tab:orange', 'tab:blue', 'tab:purple')

    data = data_dict[attribute_name]
    names = tuple(data)
    burn_times = tuple(data['EP'])
    thrusts = tuple(data['EP'][next(iter(burn_times))])
    other_names = [name for name in names if name != 'EP']

    for name in other_names:
        fig, ax = plt.subplots()
        p = []
        for burn_time, color in zip(burn_times, colors):
            for thrust, linestyle in zip(thrusts, linestyles):
                label = f'$t_b${burn_time} $F_T${thrust}'
                # Get the wanted attribute value for the current cycle
                values = data[name][burn_time][thrust]
                # Likewise for the EP-cycle and find their ratio
                ep_values = data['EP'][burn_time][thrust]
                ratios = [ep_val / val for ep_val, val in zip(ep_values, values)]
                # Plot attribute ratio and store legend info
                p_leg, = ax.plot(p_cc_range[:len(ratios)], ratios, label=label, color=color, linestyle=linestyle)
                p.append(p_leg)
        if ylim:
            ax.set_ylim(*ylim)
        base_ylabel = ylabel + r'$_{' + name.lower() + '}$ [-]'
        full_ylabel = ylabel + r'$_{ep}$ /' + base_ylabel
        ax.set_ylabel(full_ylabel)
        ax.set_xlabel('$p_{cc}$ [MPa]')
        mode = 'Frozen' if is_frozen else 'Shifting Equilibrium'
        title = 'Comparison of ' + variable_title + f' for EP-cycle and {name}-cycle'
        kwak_mode = 'KwakFix' if kwak else 'Own'
        ax.set_title(title + f'\n {mode} - {kwak_mode}')
        if line_val is not None:
            ax.plot(p_cc_range, [line_val] * len(p_cc_range), color='black', linestyle='dotted')
        if sub_legend_kwargs is None:
            sub_legend_kwargs = {}
        vertical_header_legend(burn_times, [f'{float(thrust) * 1e-3:.0f}' for thrust in thrusts], p,
                               dummy_pos=([p_cc_range[0]], [ratios[0]]),
                               **sub_legend_kwargs)
        if savefig:
            plt.savefig(dir_name + '/' + pngname + f'_{name}', dpi=1200)
        plt.show()


def value_plot(attribute_name: str, ylabel: str, variable_title: str, pngname: str, p_cc_range: tuple, is_frozen: bool,
               kwak: bool,
               data_dict: Optional[dict] = None,
               ylim: Optional[Tuple[float, float]] = None,
               sub_legend_kwargs: Optional[dict] = None,
               savefig: bool = False,
               dir_name: str = '',
               unit: str = '-',
               line_val: Optional[float] = None):
    linestyles = ('-', '--', '-.', ':')
    colors = ('tab:red', 'tab:green', 'tab:orange', 'tab:blue', 'tab:purple')

    data = data_dict[attribute_name]
    names = tuple(data)
    burn_times = tuple(data['EP'])
    thrusts = tuple(data['EP'][next(iter(burn_times))])
    for burn_time in burn_times:
        fig, ax = plt.subplots()
        p = []
        for name, color in zip(names, colors):
            for thrust, linestyle in zip(thrusts, linestyles):
                label = f'{name.lower()}-cycle $F_T${thrust}'
                # Plot attribute values and store legend info
                values = data[name][burn_time][thrust]
                p_leg, = ax.plot(p_cc_range[:len(values)], values, label=label, color=color,
                                 linestyle=linestyle)
                p.append(p_leg)
        if ylim:
            ax.set_ylim(*ylim)
        full_ylabel = ylabel + f' [{unit}]'
        ax.set_ylabel(full_ylabel)
        ax.set_xlabel('$p_{cc}$ [MPa]')
        mode = 'Frozen' if is_frozen else 'Shifting Equilibrium'
        kwak_mode = 'KwakFix' if kwak else 'Own'
        title = variable_title
        ax.set_title(title + f'\n {mode} - {kwak_mode} - $t_b$:{burn_time}')
        if sub_legend_kwargs is None:
            sub_legend_kwargs = {}
        # vertical_header_legend(names, [f'{float(thrust) * 1e-3:.0f}' for thrust in thrusts], p,
        #                        headerfunc=lambda x: x + '-cycle',
        #                        dummy_pos=([p_cc_range[0]], [values[0]]),
        #                        **sub_legend_kwargs)
        if savefig:
            plt.savefig(dir_name + '/' + pngname + f'_{burn_time}', dpi=1200)
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


def plot_dv_ratio(default_ylim=False, is_ratio: bool = True, **kwargs):
    if default_ylim:
        kwargs['ylim'] = (.72, .96)
    if is_ratio:
        plotter = ratio_plot
        ylabel = r'$I_{sp}*ln(MR^{-1})$'
        unit = 's'
    else:
        plotter = value_plot
        ylabel = r'$\Delta V$'
        unit = 'm/s'
    plotter(attribute_name='ideal_delta_v',
            ylabel=ylabel,
            variable_title=r'Ideal $\Delta$V',
            pngname='DeltaV_Comparison',
            sub_legend_kwargs={'fontsize': 7},
            unit=unit,
            **kwargs)


def plot_mr_ratio(default_ylim=False, is_ratio: bool = True, **kwargs):
    if default_ylim:
        kwargs['ylim'] = (1.2, 3.2)
    plotter = ratio_plot if is_ratio else value_plot
    plotter(attribute_name='mass_ratio',
            ylabel=r'$MR$',
            variable_title='Mass Ratio',
            pngname='Mass_Ratio_Comparison',
            sub_legend_kwargs={'loc': 2, 'fontsize': 6},
            **kwargs)


def plot_mr_kwak_ratio(default_ylim=False, is_ratio: bool = True, **kwargs):
    if default_ylim:
        kwargs['ylim'] = (1.2, 3.2)
    plotter = ratio_plot if is_ratio else value_plot
    plotter(attribute_name='mass_ratio_kwak',
            ylabel=r'$MR$',
            variable_title='Mass Ratio',
            pngname='Mass_Ratio_Comparison',
            sub_legend_kwargs={'loc': 2, 'fontsize': 6},
            **kwargs)


def plot_m0_ratio(default_ylim=False, is_ratio: bool = True, **kwargs):
    if default_ylim:
        kwargs['ylim'] = (.98, 1.010)
    plotter = ratio_plot if is_ratio else value_plot
    plotter(attribute_name='initial_mass',
            ylabel=r'm⁰',
            variable_title='Initial Mass',
            pngname='Initial_Mass_Comparison',
            line_val=1,
            unit='kg',
            **kwargs)


def plot_isp_ratio(p_cc_range: tuple, is_frozen: bool, kwak: bool,
                   data_dict: Optional[dict] = None,
                   savefig=False, dir_name: str = '',
                   default_ylim: bool = False,
                   ylim: Optional[tuple[tuple, tuple]] = None,
                   single_figure: bool = True):
    colors = {'EP': 'tab:blue',
              'GG': 'tab:red',
              'OE': 'tab:green',
              'CB': 'tab:purple'}
    ratio_colors = colors.copy() | {'GG': 'tab:cyan', 'OE': 'tab:olive', 'CB': 'tab:pink'}
    markers = {'EP': 'o',
               'GG': 's',
               'OE': '^',
               'CB': 'd'}
    ratio_markers = {'GG': 'p',
                     'OE': 'H',
                     'CB': '*'}

    data = data_dict['overall_specific_impulse']
    names = tuple(data)
    other_names = [name for name in names if name != 'EP']

    isps = {}
    for name in names:
        burn_time = next(iter(data[name]))
        thrust = next(iter(data[name][burn_time]))
        isps[name] = data[name][burn_time][thrust]

    isp_ratios = {}
    for other_name in other_names:
        isp_ratios[other_name] = [ep_isp / isp for ep_isp, isp in zip(isps['EP'], isps[other_name])]

    mode = 'Frozen' if is_frozen else 'Shifting Equilibrium'
    kwak_mode = 'KwakFix' if kwak else 'Own'
    title = 'Specific Impulse Comparison\n' + mode + ' - ' + kwak_mode
    xlabel = '$p_{cc}$ [MPa]'

    fig, ax = plt.subplots()
    ax.set_xlabel(xlabel)

    ax.set(title=title)
    if single_figure:
        ax2 = ax.twinx()
    else:
        fig2, ax2 = plt.subplots()
        ax2.set_xlabel(xlabel)
        ax2.set(title=title)
    ax.set_ylabel('$I_{sp}$ [s]')
    ax2.set_ylabel('$I_{sp}$ Ratio (EP/XX)')

    if default_ylim:
        # ax.set_ylim(340, 380)
        ax2.set_ylim(1.00, 1.08)
    elif ylim:
        ax.set_ylim(*ylim[0])
        ax2.set_ylim(*ylim[1])

    lines = []
    lines2 = []
    for name in names:
        lines += ax.plot(p_cc_range[:len(isps[name])], isps[name],
                         label=r'$I_{sp}$ ' + f'{name}-cycle',
                         color=colors[name],
                         marker=markers[name],
                         linestyle='dashed', )
    for other_name in other_names:
        lines2 += ax2.plot(p_cc_range[:len(isp_ratios[other_name])], isp_ratios[other_name],
                           label='$I_{sp}$ ratio ' + f'EP/{other_name}',
                           color=ratio_colors[other_name],
                           marker=ratio_markers[other_name],
                           markersize=7, )
    if single_figure:
        lines += lines2
    else:
        labels2 = [line.get_label() for line in lines2]
        ax2.legend(lines2, labels2)
    labels = [line.get_label() for line in lines]
    ax.legend(lines, labels)
    if savefig:
        plt.savefig(dir_name + '/' + 'Specific_Impulse_Comparison', dpi=1200)
    plt.show()


def plot_all_ratio_plots_data(short: bool = False,
                              is_frozen: bool = False,
                              burn_times: tuple = (300, 600, 900, 1200),
                              thrusts: tuple = (25e3, 50e3, 75e3, 100e3),
                              names: tuple = ('EP', 'GG'),
                              p_cc_range: tuple = tuple(range(3, 11)),
                              savefig: bool = False,
                              default_ylim: bool = False,
                              savedata: bool = False,
                              kwak: bool = False,
                              ylims: Optional[dict] = None,
                              isp_single_figure: bool = True,
                              is_ratio: bool = True,
                              **kwargs):
    if short:
        thrusts = [max(thrusts)]

    data_dict = get_data_dict(attributes=('ideal_delta_v', 'mass_ratio', 'initial_mass', 'overall_specific_impulse'),
                              thrusts=thrusts, burn_times=burn_times, names=names, kwak=kwak, is_frozen=is_frozen,
                              **kwargs)
    data_dict['info'] = kwargs | {'burn_times': burn_times, 'thrusts': thrusts, 'names': names,
                                  'is_frozen': is_frozen, 'p_cc_range': p_cc_range, 'kwak': kwak}
    if savedata:
        suffix = 'Fixed' if kwak else 'Normal'
        file_name = strftime("%Y%m%d-%H%M%S") + '_' + suffix
        with open(rf'all_plots_data\{file_name}', 'w') as f:
            json.dump(data_dict, f, indent=4)

    plot_all_ratio_from_data(data_dict=data_dict, savefig=savefig, default_ylim=default_ylim, ylims=ylims,
                             isp_single_figure=isp_single_figure, is_ratio=is_ratio)


def plot_all_ratio_from_data(data_dict: dict, savefig: bool = False, ylims: Optional[dict] = None,
                             isp_single_figure: bool = False, is_ratio: bool = True, **kwargs):
    if savefig:
        dir_name = r"all_plots_plots/" + strftime("%Y%m%d-%H%M%S")
        os.mkdir(dir_name)
        with open(dir_name + '/Matching_Data', 'w') as file:
            json.dump(data_dict, file)
    else:
        dir_name = ''

    if ylims:
        ylim_dv = ylims['dv']
        ylim_mr = ylims['mr']
        ylim_m0 = ylims['m0']
        ylim_isp = ylims['isp']
    else:
        ylim_dv, ylim_mr, ylim_m0, ylim_isp = None, None, None, None

    info = data_dict['info']
    other_kwargs = {'is_frozen': info['is_frozen'],
                    'p_cc_range': info['p_cc_range'],
                    'kwak': info['kwak'],
                    'savefig': savefig,
                    'dir_name': dir_name,
                    'data_dict': data_dict,
                    'is_ratio': is_ratio}

    plot_dv_ratio(**other_kwargs, **kwargs, ylim=ylim_dv)
    plot_mr_ratio(**other_kwargs, **kwargs, ylim=ylim_mr)
    plot_m0_ratio(**other_kwargs, **kwargs, ylim=ylim_m0)
    other_kwargs.pop('is_ratio')
    plot_isp_ratio(**other_kwargs, **kwargs, ylim=ylim_isp, single_figure=isp_single_figure)


if __name__ == '__main__':
    kwargs = {'p_cc_range': tuple(range(3, 11)),
              'is_frozen': False,
              'verbose': True,
              'kwak': False,
              'exit_pressure_forced': 0.002e6,
              'expansion_ratio_end_cooling': 20,
              'thrusts': (100e3,),
              '_ignore_cooling': True,
              'specific_impulse_quality_factor': 1,
              }

    plot_all_ratio_plots_data(default_ylim=True, savedata=True, savefig=True, names=('EP', 'GG'),
                              fuel_name='RP1_NASA', isp_single_figure=True, **kwargs)
    # def open_data_dict(filepath: str) -> dict:
    #     with open(filepath, 'r') as f:
    #         return json.load(f)
    #
    #
    # data_dict = open_data_dict(
    #     r'C:\Users\rvand\PycharmProjects\ThesisEP2\plots\KwakPlots\all_plots_data\20230108-155954_Normal')
    # plot_all_ratio_from_data(data_dict, savefig=False, isp_single_figure=False, is_ratio=False)

    # engine_data = get_engines('oe', **kwargs)
    # # print(engine_data[0].cooling_channel_section.outlet_temperature)
    # engine_data = {}
    # for name in ('EP', 'GG', 'OE'):
    #     engine_data[name] = get_engines(cycle_name=name.lower(), **kwargs)
    # plot_isp_ratio(engine_data=engine_data, p_cc_range=kwargs['p_cc_range'], is_frozen=kwargs['is_frozen'], default_ylim=False)
    # plot_all_ratio_plots(names=('EP', 'GG'), _ignore_cooling=True, verbose=True)

    # data_dict = get_data_dict(('overall_specific_impulse',), burn_times=(300,), names=('CB','EP'), **kwargs)
    # plot_isp_ratio(data_dict=data_dict, p_cc_range=kwargs['p_cc_range'], is_frozen=kwargs['is_frozen'], single_figure=True)
    # plot_all_ratio_plots_data(names=('EP', 'GG', 'OE', 'CB'), savefig=True, plot_ratios=False **kwargs)
    # with open(r'C:\Users\rvand\PycharmProjects\ThesisEP2\plots\KwakPlots\all_plots_plots\20221111-172003\Matching_Data',
    #           'r') as f:
    #     data_dict = json.load(f)
    # ylims = {'dv': (0.7, 1.0), 'mr': (1.2, 3.2), 'm0': (0.99, 1.03), 'isp': None}
    # plot_all_ratio_from_data(data_dict, savefig=True, is_ratio=True, default_ylims=True)

    # with open(r'C:\Users\rvand\PycharmProjects\ThesisEP2\plots\KwakPlots\all_plots_data\20221214-151458_Normal', 'r') as f:
    #     data_dict =json.load(f)
    # plot_all_ratio_from_data(data_dict=data_dict,savefig=False, isp_single_figure=False, is_ratio=False)

    # from plots.Imaging.mass_image import make_mass_schematic
    #
    # design_kwargs = {
    #     'is_frozen': True,
    #     'verbose': True,
    #     'exit_pressure_forced': None,
    #     'expansion_ratio': 30,
    #     'expansion_ratio_end_cooling': 20,
    #     'thrust': 100e3,
    #     'combustion_chamber_pressure': 3e6,
    #     'burn_time': 1200,
    # }
    # total_args = args.base_arguments | design_kwargs
    #
    # engine_dict = {
    #     ElectricPumpCycle: args.ep_arguments,
    #     GasGeneratorCycle: args.gg_arguments,
    #     # OpenExpanderCycle: args.oe_arguments,
    # }
    # def make_engines():
    #     for EngineClass, class_args in engine_dict.items():
    #         engine = EngineClass(**total_args, **class_args)
    #         make_mass_schematic(engine)
    #         yield engine
    # a = tuple(make_engines())
    # print()
