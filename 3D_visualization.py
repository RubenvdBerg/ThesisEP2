from EngineCycles.Abstract.EngineCycle import EngineCycle
import matplotlib.pyplot as plt
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from EngineCycles.CoolantBleedCycle import CoolantBleedCycle
from optimization import fast_optimize
from EngineArguments import arguments as args
from EngineArguments.default_arguments import get_default_kwargs
import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple, Optional
from new_new_optimization import optimize_engine

plot_defaults = {'pressure_range': (5E5, 3E7), 'mixture_range': (3.5, 6.5), 'detail_number': 10}  # LH2
plot_defaults = {'pressure_range': (5E5, 3E7), 'mixture_range': (1.5, 4.5), 'detail_number': 10}  # RP1


def set_attribute_cycle(attribute: str, cycle: EngineCycle, compare: bool = False):
    cycle_switcher = {
        ElectricPumpCycle: 'ep',
        GasGeneratorCycle: 'gg',
        CoolantBleedCycle: 'cb',
        OpenExpanderCycle: 'oe',
    }
    cycle_acronym = cycle_switcher[cycle]
    kwargs = get_default_kwargs(cycle, mass_mixture_ratio=False)

    attribute_options = ['mass', 'dv']
    if attribute == 'mass':
        symbol, unit = r'$m_0$', r'[kg]'
        title = 'Initial Mass'
    elif attribute == 'ideal_delta_v':
        symbol, unit = r'$\Delta V$', r'[m/s]'
        title = r'Ideal Change in Velocity'
    elif attribute == 'change_in_velocity_km':
        symbol, unit = r'$\Delta V$', r'[km/s]'
        title = r'Ideal Change in Velocity'
    elif attribute == 'overall_specific_impulse':
        symbol, unit = r'$I_{sp}$', '[s]'
        title = r'Specific Impulse'
    elif attribute == 'mass_ratio':
        symbol, unit = r'$\Lambda$', '[-]'
        title = 'Mass Ratio'
    else:
        raise ValueError(f'Invalid attribute. Must be one of {attribute_options}')

    normal_label = symbol + r'$_{' + cycle_acronym + '}$ ' + unit
    compare_label = symbol + r'$_{ep}$/' + normal_label

    if compare:
        return cycle_acronym, kwargs, title, compare_label
    else:
        return cycle_acronym, kwargs, title, normal_label


def threed_plot_cycle(cycle: EngineCycle,
                      pressure_range: Tuple[float, float] = plot_defaults['pressure_range'],
                      mixture_range: Tuple[float, float] = plot_defaults['mixture_range'],
                      detail_number=plot_defaults['detail_number'],
                      log: bool = True,
                      attribute='ideal_delta_v',
                      savefig: bool = False,
                      zlim: Optional[tuple] = None,
                      **kwargs):
    cycle_acronym, base_kwargs, title, z_label = set_attribute_cycle(attribute, cycle)
    kwargs = base_kwargs | kwargs
    mmr_range = np.linspace(*mixture_range, detail_number)
    pcc_range = np.logspace(*np.log10(pressure_range), detail_number) if log else np.linspace(*pressure_range,                                                                                                 detail_number)

    vmin, vmax = zlim if zlim else (None, None)
    def attributes_maker(mmr, pcc):
        try:
            print(f'MMR: {mmr:.2f}, pcc: {pcc * 1e-6:.2f} MPa')
            engine = cycle(mass_mixture_ratio=mmr,
                           combustion_chamber_pressure=pcc,
                           **kwargs)
            return getattr(engine, attribute)
        except:
            raise
            # return None

    attributes = [[attributes_maker(mmr, pcc)
                   for pcc in pcc_range]
                  for mmr in mmr_range]
    X, Y = np.meshgrid(pcc_range, mmr_range)
    Z = np.array(attributes)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    if log:
        ax.plot_surface(np.log10(X), Y, Z, cmap='viridis', edgecolor='none', vmin=vmin, vmax=vmax)
        ticks = np.log10(np.logspace(*np.log10(pressure_range), 5))
        labels = [f'{10 ** (x - 6):.1f}' for x in ticks]
        ax.set_xticks(ticks, labels)
    else:
        ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none', vmin=vmin, vmax=vmax)
        ticks = np.linspace(*pressure_range, 5)
        labels = [f'{x * 1E-6:.1f}' for x in ticks]
        ax.set_xticks(ticks, labels)
    ax.set_xlabel('$p_{cc}$ [MPa] (log)') if log else ax.set_xlabel('$p_{cc}$ [MPa]')
    ax.set_ylabel('$O/F$-ratio [-]')
    ax.set_zlabel(z_label)
    mode_name = 'Frozen at Throat' if kwargs['is_frozen'] else 'Shifting Equilibrium'
    ax.set_title(
        f'{cycle_acronym.upper()}-Cycle {title} - {mode_name} \n' + r'$F_T$=' + f'{kwargs["thrust"] * 1E-3:.1f}kN ' + r'$t_b$=' + f'{kwargs["burn_time"]:.0f}s')



    # Print min or max outputs and associated inputs
    pcc_bounds = tuple(x*1e-6 for x in pressure_range)
    x_opt, y_opt, z_opt = optimize_engine(
        CycleClass=cycle,
        total_kwargs=kwargs,
        attribute=attribute,
        bounds = (pcc_bounds,mixture_range),
    )
    z_unit = f'{z_label.split(" ")[-1].strip("[]")}'
    text = f'({x_opt:.2f} MPa, {y_opt:.2f}, {z_opt:.2f} {z_unit})'
    x_opt_val = np.log10(x_opt*1e6) if log else x_opt*1e6
    ax.set_zlim(zlim)
    ax.text(ax.get_xlim()[0], ax.get_ylim()[1], ax.get_zlim()[1]*.9, text, zdir='x', color='red')
    ax.plot(x_opt_val, y_opt, z_opt, marker='x', color='red', ls='',markersize=8, zorder=10)
    print(
        f'{cycle_acronym.upper()}-Cycle, {kwargs["thrust"] * 1e-3} kN, {kwargs["burn_time"]} s, {mode_name}: {z_opt:.4f} '
        f'{z_unit}, {x_opt:.4f} MPa, {y_opt:.4f} O/F'
    )

    if savefig:
        plt.savefig(
            'plots/3Dplots/' +
            f'Base_{attribute.upper()}_{cycle_acronym}_FT{kwargs["thrust"] * 1e-3:.0f}_tb{kwargs["burn_time"]}_{"frozen" if kwargs["is_frozen"] else "equilibrium"}'
        )
    plt.show()


def make_all_base_plots(thrust: float = 10e3, attribute: str = 'dv', **kwargs):
    for is_frozen in [True, False]:
        for cycle in ['gg', 'ep']:
            for burn_time in [300, 600, 1200]:
                threed_plot_cycle(attribute=attribute, thrust=thrust, burn_time=burn_time, cycle=cycle,
                                  is_frozen=is_frozen, **kwargs)


if __name__ == '__main__':
    from plots.KwakPlots.Results_Comparison_RP1 import engine_kwargs

    engine_kwargs.pop('combustion_chamber_pressure')
    for burn_time in [300,1200]:
        design_variables = {
            'thrust': 25e3,
            'burn_time': burn_time,
        }
        zlim = (5, 10.5) if design_variables['burn_time'] == 300 else (6.5, 11)
        total_kwargs = engine_kwargs | design_variables
        for cycle in [ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle]:
            threed_plot_cycle(cycle=cycle,
                              attribute='change_in_velocity_km',
                              log=True,
                              detail_number=3,
                              savefig=False,
                              pressure_range=(5e5, 2e7),
                              mixture_range=(1.5, 3.5),
                              zlim=zlim,
                              **total_kwargs)
