from EngineCycles.Abstract.EngineCycle import EngineCycle
import matplotlib.pyplot as plt
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle_DoublePump
from EngineCycles.CoolantBleedCycle import CoolantBleedCycle
from optimization import fast_optimize
from EngineArguments import arguments as args
from EngineArguments.default_arguments import get_default_kwargs
from plots.Imaging.performance_image import make_performance_schematic
from plots.Imaging.mass_image import make_mass_schematic
import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple

plot_defaults = {'pressure_range': (5E5, 5E6), 'mixture_range': (4, 6), 'detail_number': 100}


def set_attribute_cycle(attribute: str, cycle: EngineCycle, compare: bool = False):
    cycle_switcher = {
        ElectricPumpCycle: 'ep',
        GasGeneratorCycle: 'gg',
        CoolantBleedCycle: 'cb',
        OpenExpanderCycle_DoublePump: 'oe',
    }
    cycle_acronym = cycle_switcher[cycle]
    kwargs = get_default_kwargs(cycle)

    attribute_options = ['mass', 'dv']
    if attribute == 'mass':
        symbol, unit = r'$m_0$', r'[kg]'
        title = 'Initial Mass'
    elif attribute == 'ideal_delta_v':
        symbol, unit = r'$\Delta V$', r'[m/s]'
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


def twod_plot_cycle(cycle: EngineCycle,
                    pressure_range: Tuple[float, float] = plot_defaults['pressure_range'],
                    detail_number=plot_defaults['detail_number'],
                    attribute='ideal_delta_v',
                    savefig: bool = False,
                    best_is_min: bool = False,
                    **kwargs):
    cycle_acronym, base_kwargs, title, y_label = set_attribute_cycle(attribute, cycle)
    kwargs = base_kwargs | kwargs
    pcc_range = np.linspace(*pressure_range, detail_number)

    def attributes_maker(pcc):
        try:
            print(f'pcc: {pcc * 1e-6:.1f} MPa')
            engine = cycle(combustion_chamber_pressure=pcc,
                           **kwargs)
            # make_performance_schematic(engine)
            # make_mass_schematic(engine)
            return getattr(engine, attribute)
        except:
            return None

    attributes = [attributes_maker(pcc)
                  for pcc in pcc_range]

    fig, ax = plt.subplots()
    ax.plot(pcc_range, attributes)
    ticks = np.linspace(*pressure_range, 5)
    labels = [f'{x * 1E-6:.1f}' for x in ticks]
    ax.set_xticks(ticks, labels)
    ax.set_xlabel('$p_{cc}$ [MPa]')
    ax.set_ylabel(y_label)
    mode_name = 'Frozen at Throat' if kwargs['is_frozen'] else 'Shifting Equilibrium'
    ax.set_title(
        f'{cycle_acronym.upper()}-Cycle {title} - {mode_name} \n' + r'$F_T$=' + f'{kwargs["thrust"] * 1E-3:.1f}kN ' + r'$t_b$=' + f'{kwargs["burn_time"]:.0f}s')
    if savefig:
        plt.savefig(
            'plots/3Dplots/' +
            f'Base_{attribute.upper()}_{cycle_acronym}_FT{kwargs["thrust"] * 1e-3:.0f}_tb{kwargs["burn_time"]}_{"frozen" if kwargs["is_frozen"] else "equilibrium"}'
        )
    plt.show()

    # Print min or max outputs and associated inputs
    vals = (attribute for attribute in attributes if attribute is not None)
    best = min(vals) if best_is_min else max(vals)
    index = attributes.index(best)
    print(
        f'{cycle_acronym.upper()}-Cycle, {kwargs["thrust"] * 1e-3} kN, {kwargs["burn_time"]} s, {mode_name}: {best:.0f} '
        f'{y_label.split(" ")[-1].strip("[]")}, {pcc_range[index] * 1e-6:.4f} MPa'
    )


if __name__ == '__main__':
    design_variables = {
        'thrust': 100e3,
        'burn_time': 400,
        'verbose': False,
        'expansion_ratio': 100,
        'fuel_name': 'LH2_NASA',
        'mass_mixture_ratio': 5.5,
        'ambient_pressure': 0,
    }
    # for cycle in [ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle]:
    #     twod_plot_cycle(cycle=cycle,
    #                     attribute='overall_specific_impulse',
    #                     pressure_range=(5e5, 10e6),
    #                     **design_variables)

    twod_plot_cycle(cycle=ElectricPumpCycle,
                    attribute='mass_ratio',
                    pressure_range=(5e5, 10e6),
                    detail_number=100,
                    **design_variables)