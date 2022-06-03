from ElectricPumpCycle.EPCycle import ElectricPumpCycle
from base_gg_cycle import GasGeneratorCycle
from BaseEngineCycle.EngineCycle import EngineCycle
from typing import Optional, Tuple, Callable
import arguments as args
import matplotlib.pyplot as plt


def get_engine_plot_data(engine_cycle: type(EngineCycle), func: Callable[[type(EngineCycle)], float], thrust: float,
                         burn_time: float, p_cc_range=range(3, 11), is_frozen: bool = False, verbose: bool = False,
                         p_cc_factor=1E6, base_arguments=args.base_arguments, extra_args=args.ep_arguments):
    engines = [engine_cycle(combustion_chamber_pressure=p_cc * p_cc_factor,
                            thrust=thrust,
                            burn_time=burn_time,
                            is_frozen=is_frozen,
                            verbose=verbose,
                            **base_arguments, **extra_args) for p_cc in p_cc_range]
    return [func(engine) for engine in engines]


def plot_isp_ratio(thrust: float, burn_time: float, is_frozen: bool = False, verbose: bool = False, savefig=False):
    p_ccs = range(3, 11)
    gg_isps = get_engine_plot_data(GasGeneratorCycle, lambda x: x.simple_specific_impulse, extra_args=args.gg_arguments)
    ep_isps = get_engine_plot_data(ElectricPumpCycle, lambda x: x.simple_specific_impulse, extra_args=args.ep_arguments)
    isp_ratios = [ep / gg for ep, gg in zip(ep_isps, gg_isps)]
    fig, ax = plt.subplots()
    ax.set_xlabel('$p_{cc}$ [MPa]')
    ax.set_ylabel('$I_{sp}$ [s]')
    ax.set(title='Specific impulse comparison')
    ep_line = ax.plot(p_ccs, gg_isps, label='$I_{sp}$ (GG Cycle)', color='tab:red', linestyle='dashed')
    gg_line = ax.plot(p_ccs, ep_isps, label='$I_{sp}$ (EP Cycle)', color='tab:blue', linestyle='dashed')
    ax.set_ylim(340, 380)
    ax2 = ax.twinx()
    ax2.set_ylabel('$I_{sp}$ Ratio (EP/GG)')
    ratio_line = ax2.plot(p_ccs, isp_ratios, label='$I_{sp}$ Ratio', color='#29c6cf')
    ax2.set_ylim(1.00, 1.08)
    lines = gg_line + ep_line + ratio_line
    labels = [line.get_label() for line in lines]
    ax.legend(lines, labels, loc=2)
    if savefig:
        plt.savefig('Specific_Impulse_Comparison.png')
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


def ratio_plot(func, ylabel, title, pngname, engines_dict={},
               plot_straight_line_at=0,
               p_cc_range=range(3, 11),
               ylim: Optional[Tuple[float, float]] = None,
               xlabel='$p_{cc}$ [MPa]',
               burn_times=(300, 600, 900, 1200),
               thrusts=(25, 50, 75, 100),
               linestyles=('-', '--', '-.', ':'),
               colors=('tab:red', 'tab:green', 'tab:orange', 'tab:blue', 'tab:purple'),
               short=False,
               sub_legend_kwargs={},
               savefig=False,
               ):
    if short:
        thrusts = [thrusts[-1]]
    p = []
    fig, ax = plt.subplots()
    for burn_time, colour in zip(burn_times, colors):
        for thrust, linestyle in zip(thrusts, linestyles):
            arguments['thrust'] = thrust * 1E3
            arguments['burn_time'] = burn_time
            if engines_dict:
                ep_datas = [func(engine) for engine in engines_dict[f'EP tb:{burn_time}, ft:{thrust}']]
                gg_datas = [func(engine) for engine in engines_dict[f'GG tb:{burn_time}, ft:{thrust}']]
            else:
                ep_datas = get_engine_plot_data(ElectricPumpCycle, func, args=arguments, extra_args=ep_arguments)
                gg_datas = get_engine_plot_data(GasGeneratorCycle, func, args=arguments, extra_args=gg_arguments)
            data_ratios = [ep_data / gg_data for ep_data, gg_data in zip(ep_datas, gg_datas)]
            p_leg, = ax.plot(p_cc_range, data_ratios,
                             label=f'$t_b${burn_time} $F_T${burn_time}',
                             color=colour,
                             linestyle=linestyle)
            p.append(p_leg)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    mode = 'Frozen' if arguments['is_frozen'] else 'Shifting Equilibrium'
    ax.set_title(title + f'\n {mode}')
    if plot_straight_line_at:
        line_val = plot_straight_line_at
        ax.plot(p_cc_range, [line_val] * len(p_cc_range), color='black', linestyle='dotted')
    vertical_header_legend(burn_times, thrusts, p, **sub_legend_kwargs)
    if savefig:
        plt.savefig(pngname)
    plt.show()


def plot_dv_ratio(default_vals=True, **kwargs):
    if default_vals:
        kwargs['ylim'] = (.72, .96)
    dv_title = 'Comparison of ideal DeltaV of EP Cycle to GG Cycle'
    dv_ylabel = '$(I_{sp}*ln(MR^{-1}))_{ep}$ / $(I_{sp}*ln(MR^{-1}))_{gg}$ [-]'
    dv_legend_kwargs = {'fontsize': 7}
    dv_pngname = 'DeltaV_Comparison'
    ratio_plot(lambda x: x.ideal_delta_v, dv_ylabel, dv_title, dv_pngname,
               sub_legend_kwargs=dv_legend_kwargs, **kwargs)


def plot_mr_ratio(default_vals=True, **kwargs):
    if default_vals:
        kwargs['ylim'] = (1.2, 3.2)
    mr_title = 'Comparison of EP Cycle mass ratio to GG Cycle mass ratio'
    mr_ylabel = '$MR_{ep}/MR_{gg}$ [-]'
    mr_legend_kwargs = {'dummy_pos': ([3], [1.5]), 'loc': 2, 'fontsize': 6}
    mr_pngname = 'Mass_Ratio_Comparison'
    ratio_plot(lambda x: x.mass_ratio, mr_ylabel, mr_title, mr_pngname,
               sub_legend_kwargs=mr_legend_kwargs, **kwargs)


def plot_m0_ratio(default_vals=True, **kwargs):
    if default_vals:
        kwargs['ylim'] = (.98, 1.010)
    m0_title = 'Initial mass comparison: EP Cycle to GG Cycle'
    m0_ylabel = r'$m⁰_{ep}/m⁰_{gg}$ [-]'
    m0_pngname = 'Initial_Mass_Comparison'
    ratio_plot(lambda x: x.mass, m0_ylabel, m0_title, m0_pngname, plot_straight_line_at=1, **kwargs)


def plot_all_ratio_plots(short: bool = False, **kwargs):
    kwargs['short'] = short
    burn_times = (300, 600, 900, 1200)
    thrusts = [100] if short else (25, 50, 75, 100)
    p_cc_range = range(3, 11)
    labels = ('EP', 'GG')
    classes = (ElectricPumpCycle, GasGeneratorCycle)
    extra_kwargses = (ep_arguments, gg_arguments)
    enginedict = {}
    for label, Class, extra_kwargs in zip(labels, classes, extra_kwargses):
        for burn_time in burn_times:
            for thrust in thrusts:
                arguments['burn_time'] = burn_time
                arguments['thrust'] = thrust
                enginedict[f'{label} tb:{burn_time}, ft:{thrust}'] = [
                    Class(combustion_chamber_pressure=p_cc * 1E6, **arguments, **extra_kwargs)
                    for p_cc in p_cc_range
                ]
    plot_isp_ratio()
    plot_dv_ratio(engines_dict=enginedict, **kwargs)
    plot_mr_ratio(engines_dict=enginedict, **kwargs)
    plot_m0_ratio(engines_dict=enginedict, **kwargs)


if __name__ == '__main__':
    kwargs = {'savefig': True,
              'default_vals': True}
    # plot_isp_ratio()
    # plot_mr_ratio(**kwargs)
    # plot_dv_ratio(**kwargs)
    # plot_m0_ratio(**kwargs)
    plot_all_ratio_plots(short=True)
