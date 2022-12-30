from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from optimization import fast_optimize
from EngineArguments import arguments as args
import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple

base_args = args.base_arguments_o
plot_defaults = {'pressure_range': (5E5, 3E7), 'mixture_range': (1.5, 4.0), 'detail_number': 50}


def set_attribute_cycle(attribute: str, cycle_type: str, compare: bool = False):
    cycle_types = ['ep', 'gg', 'ex', 'sc']
    if cycle_type == 'ep':
        cycle = ElectricPumpCycle
        extra_arguments = args.ep_arguments
    elif cycle_type == 'ex':
        raise NotImplementedError
    elif cycle_type == 'gg':
        cycle = GasGeneratorCycle
        extra_arguments = args.gg_arguments
    elif cycle_type == 'sc':
        raise NotImplementedError
    else:
        raise ValueError(f'Invalid cycle_type. Must be one of {cycle_types}')

    attribute_options = ['mass', 'dv']
    if attribute == 'mass':
        def attribute_function(x):
            return x.mass

        z_label = r'$m_{0_{' + cycle_type + r'}}$ [kg]'
        z_label_compare = r'$m_{0_{ep}}$/$m_{0_{' + cycle_type + r'}}$ [kg]'
        title = 'Initial Mass'
    elif attribute == 'dv':
        def attribute_function(x):
            return x.ideal_delta_v

        z_label = r'$\Delta V_{' + cycle_type + r'}$ [km/s]'
        z_label_compare = r'$\Delta V_{ep}$/$\Delta V_{' + cycle_type + r'}$ [km/s]'
        title = r'Ideal $\Delta$V'
    else:
        raise ValueError(f'Invalid attribute. Must be one of {attribute_options}')

    if compare:
        return cycle, extra_arguments, title, attribute_function, z_label_compare
    else:
        return cycle, extra_arguments, title, attribute_function, z_label


def threed_plot_cycle(thrust: float, burn_time: float, is_frozen: bool = False, verbose: bool = False,
                      pressure_range: Tuple[float, float] = plot_defaults['pressure_range'], log: bool = True,
                      mixture_range: Tuple[float, float] = plot_defaults['mixture_range'], cycle_type='ep',
                      detail_number=plot_defaults['detail_number'],
                      attribute='dv', base_arguments: dict = base_args, savefig: bool = False):
    cycle, extra_arguments, title, attribute_function, z_label = set_attribute_cycle(attribute, cycle_type)
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
    if savefig:
        plt.savefig(
            'plots/3Dplots/' +
            f'Base_{attribute.upper()}_{cycle_type}_FT{thrust * 1e-3:.0f}_tb{burn_time}_{"frozen" if is_frozen else "equilibrium"}'
        )
    plt.show()

    # Print min or max outputs and associated inputs
    a, arg = (np.amax, np.argmax) if attribute == 'dv' else (np.amin, np.argmin)
    Z_ult = a(Z)
    Z_index = np.unravel_index(arg(Z), np.shape(Z))
    print(
        f'{cycle_type.upper()}-Cycle, {thrust * 1e-3} kN, {burn_time} s, {mode_name}: {Z_ult:.0f} '
        f'{z_label.split(" ")[-1].strip("[]")}, {X[Z_index] * 1e-6:.4f} MPa, {Y[Z_index]:.4f} O/F'
    )


def threed_plot_comparison_cycle(thrust: float, burn_time: float, is_frozen: bool = False, verbose: bool = False,
                                 pressure_range: Tuple[float, float] = plot_defaults['pressure_range'],
                                 log: bool = True, mixture_range: Tuple[float, float] = plot_defaults['mixture_range'],
                                 cycle_type='gg', detail_number=20, attribute='mass', base_arguments: dict = base_args,
                                 savefig: bool = False):
    cycle, extra_arguments, title, attribute_function, z_label = set_attribute_cycle(attribute, cycle_type,
                                                                                     compare=True)
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
                                                           thrust=thrust,
                                                           burn_time=burn_time,
                                                           is_frozen=is_frozen,
                                                           verbose=verbose,
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
    if savefig:
        plt.savefig(
            'plots/3Dplots/' +
            fr'Comparison_{attribute.upper()}_ep_{cycle_type}_FT{thrust * 1e-3:.0f}_tb{burn_time}_{"frozen" if is_frozen else "equilibrium"}_{"log" if log else ""}'
        )
    plt.show()


def threed_plot_cycle_opt(is_frozen: bool = False, verbose: bool = False,
                          thrust_range: Tuple[float, float] = (1e3, 100e3),
                          burn_time_range: Tuple[float, float] = (300, 1200), log: bool = False, cycle_type='ep',
                          detail_number=25, attribute='dv', savefig: bool = False):
    cycle, extra_arguments, title, attribute_function, z_label = set_attribute_cycle(attribute, cycle_type)

    tb_range = np.linspace(*burn_time_range, detail_number)
    ft_range = np.logspace(*np.log10(thrust_range), detail_number) if log else np.linspace(*thrust_range,
                                                                                           detail_number)
    opt_variables = [[fast_optimize(thrust=thrust,
                                    burn_time=burn_time,
                                    is_frozen=is_frozen,
                                    verbose=verbose,
                                    cycle_type=cycle_type)
                      for thrust in ft_range]
                     for burn_time in tb_range]
    X, Y = np.meshgrid(ft_range, tb_range)
    Z = np.array(opt_variables)
    if attribute == 'dv':
        Z = Z / 1e3

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    if log:
        ax.plot_surface(np.log10(X), Y, Z, cmap='viridis', edgecolor='none')
        ticks = np.log10(np.logspace(*np.log10(thrust_range), 5))
        labels = [f'{10 ** (x - 3):.1f}' for x in ticks]
        ax.set_xticks(ticks, labels)
    else:
        ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')
        ticks = np.linspace(*thrust_range, 5)
        labels = [f'{x * 1E-3:.1f}' for x in ticks]
        ax.set_xticks(ticks, labels)
    ax.set_xlabel('$F_{T}$ [kN] (log)') if log else ax.set_xlabel('$F_{T}$ [kN]')
    ax.set_ylabel('$t_b$ [s]')
    ax.set_zlabel(z_label, labelpad=8)
    mode_name = 'Frozen at Throat' if is_frozen else 'Shifting Equilibrium'
    ax.set_title(f'{cycle_type.upper()}-Cycle Optimized {title} - {mode_name}')
    plt.show()
    if savefig:
        plt.savefig(
            'plots/3Dplots/' +
            fr'Base_Optimized_{attribute.upper()}__{cycle_type}_{"frozen" if is_frozen else "equilibrium"}_{"log" if log else ""}'
        )


def make_all_base_plots(thrust: float = 10e3, attribute: str = 'dv', **kwargs):
    for is_frozen in [True, False]:
        for cycle in ['gg', 'ep']:
            for burn_time in [300, 600, 1200]:
                threed_plot_cycle(attribute=attribute, thrust=thrust, burn_time=burn_time, cycle_type=cycle,
                                  is_frozen=is_frozen, **kwargs)


def make_all_comparison_plots(thrust: float = 10e3, attribute: str = 'dv', **kwargs):
    for is_frozen in [True, False]:
        for burn_time in [300, 600, 1200]:
            threed_plot_comparison_cycle(attribute=attribute, thrust=thrust, burn_time=burn_time, is_frozen=is_frozen,
                                         **kwargs)


if __name__ == '__main__':
    # make_all_base_plots(savefig=True, attribute='mass')
    # make_all_comparison_plots(savefig=True, attribute='mass', log=False)
    # threed_plot_cycle_opt(cycle_type='ep',
    #                       attribute='dv',
    #                       log=False,
    #                       detail_number=20,
    #                       savefig=True,
    #                       is_frozen=False)
    threed_plot_cycle(10e3, 300)
    # threed_plot_comparison_cycle(thrust=10e3, burn_time=300, cycle_type='gg', savefig=True)

    # threed_plot_cycle(thrust=10E3,
    #                   burn_time=300,
    #                   pressure_range=(1E5, 3E7),
    #                   mixture_range=(1.5, 4.0),
    #                   cycle_type='ep',
    #                   attribute='dv',
    #                   detail_number=50,
    #                   log=True,
    #                   is_frozen=False)
    # threed_plot_cycle_opt(detail_number=2, attribute='dv')
