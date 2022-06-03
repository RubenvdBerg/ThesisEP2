from ElectricPumpCycle.EPCycle import ElectricPumpCycle
import arguments as args
import matplotlib.pyplot as plt

arguments = {'exit_pressure': .002E6, 'max_acceleration': 4.5, 'heat_ratio_pressurant': 1.667,
             'mass_mixture_ratio': 2.45, 'pressurant_initial_pressure': 27E6, 'pressurant_final_pressure': 5E6,
             'oxidizer_initial_pressure': .4E6, 'fuel_initial_pressure': .25E6, 'fuel_pump_pressure_factor': 1.55,
             'oxidizer_pump_pressure_factor': 1.15, 'pressurant_gas_constant': 2080,
             'pressurant_initial_temperature': 100, 'oxidizer_pump_efficiency': .66, 'fuel_pump_efficiency': .61,
             'fuel_pump_specific_power': 15E3, 'oxidizer_pump_specific_power': 20E3, 'pressurant_margin_factor': 1.1,
             'pressurant_tank_structural_factor': 1.2, 'propellant_margin_factor': 1.01, 'tanks_structural_factor': 2.5,
             'ullage_volume_factor': 1.08, 'oxidizer_density': 1126.1, 'fuel_density': 804.2,
             'tanks_material_density': 2850, 'pressurant_tank_material_density': 4430, 'tanks_yield_strength': 250E6,
             'pressurant_tank_yield_strength': 1100E6, 'fuel_specific_heat': 2009,
             'electric_motor_specific_power': 5.3E3, 'inverter_specific_power': 60E3, 'battery_specific_power': 6.95E3,
             'battery_specific_energy': 198 * 3600, 'electric_motor_efficiency': .95, 'inverter_efficiency': .85,
             'battery_structural_factor': 1.2, 'coolant_allowable_temperature_change': 40, 'verbose': False}

arguments = args.base_arguments | args.ep_arguments | {'is_frozen': False}

burn_times_def = (300, 390, 1200)
# ylims = [[(7.55, 9.5), (0, 1.6)], [(7.55, 9.5), (0, 1.6)], [(7.55, 9.5), (0, 1.6)]]  # Consistent axes
ylims_def = ([(7.2, 9.5), (0, 2.3)], [(7.55, 9.5), (0, 1.6)], [(7.55, 8.6), (0, .099)])  # Axes equal to paper
chamber_pressures_def = range(3, 11)


def broken_axis(plotfunc, y_limits):
    figre, (axis_1, axis_2) = plt.subplots(2, 1)
    plotfunc(axis_1)
    plotfunc(axis_2)
    axis_1.set_ylim(*y_limits[0])
    axis_2.set_ylim(*y_limits[1])
    axis_1.spines.bottom.set_visible(False)
    axis_2.spines.top.set_visible(False)
    axis_1.xaxis.tick_top()
    axis_1.tick_params(labeltop=False)  # Don't put tick labels at the top
    axis_2.xaxis.tick_bottom()
    d = .5
    kwargs = dict(
        marker=[(-1, -d), (1, d)], markersize=9, linestyle='none', color='k', mec='k', mew=1, clip_on=False
    )
    axis_1.plot([0, 1], [0, 0], transform=axis_1.transAxes, **kwargs)
    axis_2.plot([0, 1], [1, 1], transform=axis_2.transAxes, **kwargs)
    return figre, (axis_1, axis_2)


def plot_fuel_flow_battery(burn_times=burn_times_def, ylimits=ylims_def, chamber_pressures=chamber_pressures_def,
                           savefig=False):
    for burn_time, ylim in zip(burn_times, ylimits):
        print(f'Making engines for tb:{burn_time:.0f}')
        engines = [ElectricPumpCycle(
            thrust=100E3,
            burn_time=burn_time,
            combustion_chamber_pressure=p_cc * 1E6,
            **arguments
        )
            for p_cc in chamber_pressures]
        total_fuel_flow = [engine.total_fuelpump_flow for engine in engines]
        coolant_flow = [engine.battery.coolant_flow_required for engine in engines]
        main_flow = [engine.total_fuelpump_flow - engine.battery.coolant_flow_required for engine in engines]
        for main, checkengine in zip(main_flow, engines):
            if main != checkengine.fuel.mass_flow:
                print(f'main and check should be equal, but they are: {main}, {checkengine.fuel.mass_flow} '
                      f'for t_b {burn_time} s, pcc {checkengine.combustion_chamber_pressure / 1E6} MPa')

        def plots(axiss):
            axiss.plot(chamber_pressures, total_fuel_flow, label='Fuel Pump Mass Flow w/ Regenerative Cooling')
            axiss.plot(chamber_pressures, main_flow, label='Fuel Pump Mass Flow w/o Regenerative Cooling')
            axiss.plot(chamber_pressures, coolant_flow, label='Coolant Mass Flow')

        fig, (ax1, ax2) = broken_axis(plots, ylim)
        fig.text(.04, .5, 'Fuel Mass Flow Rate [kg/s]', va='center', rotation='vertical')
        ax2.set_xlabel('$p_{cc}$ [MPaA]')
        ax1.set(title=f'Fuel mass flow rate for $F_T=100kN$ and $t_b={burn_time} s$')
        plt.legend(loc=3, bbox_to_anchor=(0, 0.85))
        if savefig:
            plt.savefig(f'Fuel_Mass_Flow_EP_{burn_time}.png')
        plt.show()


if __name__ == '__main__':
    plot_fuel_flow_battery()
