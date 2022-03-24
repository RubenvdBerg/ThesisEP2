from base_ep_cycle import ElectricPumpCycle
from base_gg_cycle import GasGeneratorCycle
from optimization import InitialMassOpt, DeltaVOpt, optimal_cycle_variables
import arguments as args
import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple

base_args = args.base_arguments
del base_args['mass_mixture_ratio']


def threed_plot_cycle(thrust: float, burn_time: float, is_frozen: bool = False, verbose: bool = False,
                      pressure_range: Tuple[float, float] = (1E6, 20E6), log: bool = True,
                      mixture_range: Tuple[float, float] = (1.8, 3.5), cycle_type='ep', detail_number=20,
                      attribute='mass', base_arguments: dict = base_args):
    cycle_types = ['ep', 'gg', 'ex', 'sc']
    assert cycle_type in cycle_types, f'Invalid cycle_type. Must be one of {cycle_types}'
    if cycle_type == 'ep':
        extra_arguments = args.ep_arguments
        cycle = ElectricPumpCycle
    elif cycle_type == 'ex':
        raise NotImplementedError
    elif cycle_type == 'gg':
        extra_arguments = args.gg_arguments
        cycle = GasGeneratorCycle
    elif cycle_type == 'sc':
        raise NotImplementedError

    attribute_options = ['mass', 'dv']
    assert attribute in attribute_options, f'Invalid attribute. Must be one of {attribute_options}'
    if attribute == 'mass':
        def attribute_function(x): return x.mass

        z_label = r'$m_{0_{' + cycle_type + r'}}$ [kg]'
        title = 'Initial Mass'
    if attribute == 'dv':
        def attribute_function(x): return x.ideal_delta_v

        z_label = r'$\Delta V_{' + cycle_type + r'}$ [m/s]'
        title = r'Ideal $\Delta$V'

    mmr_range = np.linspace(*mixture_range, detail_number)
    pcc_range = np.logspace(*np.log10(pressure_range), detail_number) if log else np.linspace(*pressure_range,
                                                                                              detail_number)

    attributes = [[attribute_function(cycle(mass_mixture_ratio=mmr,
                                            combustion_chamber_pressure=pcc,
                                            thrust=thrust,
                                            burn_time=burn_time,
                                            is_frozen=is_frozen,
                                            verbose=verbose,
                                            **base_arguments,
                                            **extra_arguments
                                            ))
                   for pcc in pcc_range]
                  for mmr in mmr_range]
    X, Y = np.meshgrid(pcc_range, mmr_range)
    Z = np.array(attributes)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    if log:
        ax.plot_surface(np.log10(X), Y, Z, cmap='viridis', edgecolor='none')
        ticks = np.log10(np.logspace(*np.log10(pressure_range), 5))
        labels = [f'{10 ** (x - 6):.1f}' for x in ticks]
        ax.set_xticks(ticks, labels)
    else:
        ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')
        ticks = np.linspace(*pressure_range, 5)
        labels = [f'{x * 1E-6:.1f}' for x in ticks]
        ax.set_xticks(ticks, labels)
    ax.set_xlabel('$p_{cc}$ [MPa] (log)') if log else ax.set_xlabel('$p_{cc}$ [MPa]')
    ax.set_ylabel('$O/F$-ratio [-]')
    ax.set_zlabel(z_label)
    mode_name = 'Frozen at Throat' if is_frozen else 'Shifting Equilibrium'
    ax.set_title(
        f'{cycle_type.upper()}-Cycle {title} - {mode_name} \n' + r'$F_T$=' + f'{thrust * 1E-3:.1f}kN ' + r'$t_b$=' + f'{burn_time:.0f}s')
    plt.show()


def threed_plot_comparison_cycle(thrust: float, burn_time: float, is_frozen: bool = False, verbose: bool = False,
                                 pressure_range: Tuple[float, float] = (1E6, 20E6), log: bool = True,
                                 mixture_range: Tuple[float, float] = (1.8, 3.5), cycle_type='ep', detail_number=20,
                                 attribute='mass', base_arguments: dict = base_args):
    cycle_types = ['gg', 'ex', 'sc']
    assert cycle_type in cycle_types, f'Invalid cycle_type. Must be one of {cycle_types}'
    if cycle_type == 'ex':
        raise NotImplementedError
    elif cycle_type == 'gg':
        extra_arguments = args.gg_arguments
        cycle = GasGeneratorCycle
    elif cycle_type == 'sc':
        raise NotImplementedError

    attribute_options = ['mass', 'dv']
    assert attribute in attribute_options, f'Invalid attribute. Must be one of {attribute_options}'
    if attribute == 'mass':
        def attribute_function(x): return x.mass

        z_label = r'$m_{0_{ep}}$/$m_{0_{' + cycle_type + r'}}$ [kg]'
        title = 'Initial Mass'
    if attribute == 'dv':
        def attribute_function(x): return x.ideal_delta_v

        z_label = r'$\Delta V_{ep}$/$\Delta V_{' + cycle_type + r'}$ [m/s]'
        title = r'Ideal $\Delta$V'

    mmr_range = np.linspace(*mixture_range, detail_number)
    pcc_range = np.logspace(*np.log10(pressure_range), detail_number) if log else np.linspace(*pressure_range,
                                                                                              detail_number)

    compare_attributes = [[attribute_function(cycle(mass_mixture_ratio=mmr,
                                                    combustion_chamber_pressure=pcc,
                                                    thrust=thrust,
                                                    burn_time=burn_time,
                                                    is_frozen=is_frozen,
                                                    verbose=verbose,
                                                    **base_arguments,
                                                    **extra_arguments
                                                    ))
                           for pcc in pcc_range]
                          for mmr in mmr_range]
    ep_attributes = [[attribute_function(ElectricPumpCycle(mass_mixture_ratio=mmr,
                                                           combustion_chamber_pressure=pcc,
                                                           **base_arguments,
                                                           **args.ep_arguments
                                                           ))
                      for pcc in pcc_range]
                     for mmr in mmr_range]
    attribute_ratios = np.array(ep_attributes) / np.array(compare_attributes)

    X, Y = np.meshgrid(pcc_range, mmr_range)
    Z = np.array(attribute_ratios)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    if log:
        ax.plot_surface(np.log10(X), Y, Z, cmap='viridis', edgecolor='none')
        ticks = np.log10(np.logspace(*np.log10(pressure_range), 5))
        labels = [f'{10 ** (x - 6):.1f}' for x in ticks]
        ax.set_xticks(ticks, labels)
    else:
        ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')
        ticks = np.linspace(*pressure_range, 5)
        labels = [f'{x * 1E-6:.1f}' for x in ticks]
        ax.set_xticks(ticks, labels)
    ax.set_xlabel('$p_{cc}$ [MPa] (log)') if log else ax.set_xlabel('$p_{cc}$ [MPa]')
    ax.set_ylabel('$O/F$-ratio [-]')
    ax.set_zlabel(z_label)
    mode_name = 'Frozen at Throat' if is_frozen else 'Shifting Equilibrium'
    ax.set_title(
        f'Comparison {cycle_type.upper()}- and EP-Cycle {title} - {mode_name} \n' + r'$F_T$=' + f'{thrust * 1E-3:.1f}kN ' + r'$t_b$=' + f'{burn_time:.0f}s')
    plt.show()


def threed_plot_cycle_opt(opt_method: str = 'Nelder-Mead', n_pop: int = 1, is_frozen: bool = False,
                          verbose: bool = False, thrust_range: Tuple[float, float] = (1e3, 100e3),
                          burn_time_range: Tuple[float, float] =(300, 1200), log: bool = False, cycle_type='ep',
                          detail_number=20, attribute='mass'):
    cycle_types = ['ep', 'gg', 'ex', 'sc']
    assert cycle_type in cycle_types, f'Invalid cycle_type. Must be one of {cycle_types}'
    if cycle_type == 'ep':
        cycle = ElectricPumpCycle
    elif cycle_type == 'ex':
        raise NotImplementedError
    elif cycle_type == 'gg':
        cycle = GasGeneratorCycle
    elif cycle_type == 'sc':
        raise NotImplementedError

    attribute_options = ['mass', 'dv']
    assert attribute in attribute_options, f'Invalid attribute. Must be one of {attribute_options}'
    if attribute == 'mass':
        opt_class = InitialMassOpt
        z_label = r'$m_{0_{' + cycle_type + r'}}$ [kg]'
        title = 'Initial Mass'
    if attribute == 'dv':
        opt_class = DeltaVOpt
        z_label = r'$\Delta V_{' + cycle_type + r'}$ [m/s]'
        title = r'Ideal $\Delta$V'

    tb_range = np.linspace(*burn_time_range, detail_number)
    thrust_range = np.logspace(*np.log10(thrust_range), detail_number) if log else np.linspace(*thrust_range,
                                                                                               detail_number)
    opt_variables = [[optimal_cycle_variables(thrust=thrust,
                                                                                         burn_time=burn_time,
                                                                                         is_frozen=is_frozen,
                                                                                         verbose=verbose,
                                                                                         cycle=cycle,
                                                                                         method=opt_method,
                                                                                         n_pop=n_pop,
                                                                                         optimize_class=opt_class)
                   for thrust in thrust_range]
                  for burn_time in tb_range]
    attributes = [opt_outputs[0] for opt_inputs, opt_outputs in opt_variables]
    X, Y = np.meshgrid(thrust_range, tb_range)
    Z = np.array(attributes)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    if log:
        ax.plot_surface(np.log10(X), Y, Z, cmap='viridis', edgecolor='none')
        ticks = np.log10(np.logspace(*np.log10(thrust_range), 5))
        labels = [f'{10 ** (x - 6):.1f}' for x in ticks]
        ax.set_xticks(ticks, labels)
    else:
        ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')
        ticks = np.linspace(*thrust_range, 5)
        labels = [f'{x * 1E-6:.1f}' for x in ticks]
        ax.set_xticks(ticks, labels)
    ax.set_xlabel('$F_{T}$ [N] (log)') if log else ax.set_xlabel('$F_{T}$ [N]')
    ax.set_ylabel('$t_b$ [s]')
    ax.set_zlabel(z_label)
    mode_name = 'Frozen at Throat' if is_frozen else 'Shifting Equilibrium'
    ax.set_title(f'{cycle_type.upper()}-Cycle Optimized {title} - {mode_name}')
    plt.show()


# "2.71271150e+00, 1.36836975e+06]), array([-11642.43734275])" 1200
# "2.65939854e+00, 1.01274482e+06]), array([-11521.10834996]" 300

if __name__ == '__main__':
    # for cycle_type in ['gg', 'ep']:
    #     for burn_time in [300, 1200]:
    #         base_arguments['burn_time'] = burn_time
    #         threed_plot_cycle(cycle_type=cycle_type, arguments=base_arguments, detail_number=15, log=True,
    #                           attribute=lambda x: x.mass_ratio, title='Mass Ratio')
    # threed_plot_comparison_cycle(cycle_type='gg', attribute=lambda x: x.mass, title='Initial Mass', detail_number=5, log=False)
    # threed_plot_comparison_cycle(cycle_type='gg', attribute='mass', detail_number=5, log=False)
    threed_plot_cycle(thrust=10E3,
                      burn_time=300,
                      pressure_range=(1E5, 20E6),
                      cycle_type='ep',
                      attribute='dv',
                      detail_number=15,
                      log=True)
    # threed_plot_cycle_opt(detail_number=2, attribute='dv')