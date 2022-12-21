import datetime

import pandas as pd
import arguments as args
from scipy import constants
from EngineCycles.GasGeneratorCycle.GGCycle import GasGeneratorCycle
from EngineCycles.ElectricPumpCycle.EPCycle import ElectricPumpCycle
from KwakFix.KwakFixCycles import KwakFixElectricPumpCycle, KwakFixGasGeneratorCycle
from itertools import zip_longest

burn_times = (300, 390, 1200)


# batt_coolant_flows = (0.979789955909310, 0.853084378343365, 0.0738446782915450)


def full_output(kwak_fix, design_args=args.desgin_arguments, common_args=args.common_arguments_kwak,):
    combined_args = design_args | common_args | {'exit_pressure_forced': 0.002E6}
    if kwak_fix:
        combined_args = combined_args | {'fuel_pump_pressure_factor': 1.55, 'oxidizer_pump_pressure_factor': 1.15, }
    if kwak_fix:
        ep_args = args.ep_arguments_rp1_kwak
        gg_args = args.gg_arguments_rp1_kwak | {'_exhaust_thrust_contribution': .01}
        gg_cycle_class = KwakFixGasGeneratorCycle
        ep_cycle_class = KwakFixElectricPumpCycle
    else:
        ep_args = args.ep_arguments
        gg_args = args.gg_arguments
        gg_cycle_class = GasGeneratorCycle
        ep_cycle_class = ElectricPumpCycle
    gg_cycle = gg_cycle_class(**combined_args, **gg_args, iterate=False)
    ep_cycle = ep_cycle_class(**combined_args, **ep_args, iterate=False)
    arguments = args.common_arguments_kwak
    input_col1 = [
        ep_cycle.battery_cooler.coolant_specific_heat_capacity,
        gg_cycle.gas_generator.outlet_flow_state.specific_heat_capacity,
        arguments['max_acceleration'],
        arguments['pressurant_heat_capacity_ratio'],
        gg_cycle.gas_generator.outlet_flow_state.heat_capacity_ratio,
        arguments['mass_mixture_ratio'],
        gg_cycle.gas_generator.mass_mixture_ratio,
        arguments['pressurant_initial_pressure'] * 1E-6,
        arguments['pressurant_final_pressure'] * 1E-6,
        arguments['fuel_initial_pressure'] * 1E-6,
        arguments['oxidizer_initial_pressure'] * 1E-6,
        gg_args['turbine_pressure_ratio'],
        constants.R / arguments['pressurant_molar_mass'],
        gg_cycle.gas_generator.outlet_flow_state.specific_gas_constant,
        arguments['pressurant_initial_temperature'],
        gg_args['turbine_maximum_temperature'],
        gg_args['gg_stay_time'] * 1E3,
        gg_args['turbopump_specific_power'] * 1E-3,
        ep_args['oxidizer_pump_specific_power'] * 1E-3,
        ep_args['fuel_pump_specific_power'] * 1E-3,
        ep_args['electric_motor_specific_power'] * 1E-3,
        ep_args['inverter_specific_power'] * 1E-3,
        ep_args['battery_specific_power'] * 1E-3
    ]
    input_col2 = [
        ep_args['battery_specific_energy'] / 3600,
        arguments['oxidizer_pump_efficiency'],
        arguments['fuel_pump_efficiency'],
        gg_args['turbine_efficiency'],
        ep_args['electric_motor_efficiency'],
        ep_args['inverter_efficiency'],
        ep_args['battery_structural_factor'],
        arguments['pressurant_margin_factor'],
        gg_args['gg_structural_factor'],
        arguments['pressurant_tank_safety_factor'],
        arguments['propellant_margin_factor'],
        arguments['tanks_structural_factor'],
        arguments['ullage_volume_factor'],
        arguments['oxidizer_density'],
        arguments['fuel_density'],
        gg_args['gg_material_density'],
        arguments['tanks_material_density'],
        arguments['tanks_material_density'],
        arguments['pressurant_tank_material_density'],
        gg_args['gg_yield_strength'] * 1E-6,
        arguments['tanks_yield_strength'] * 1E-6,
        arguments['tanks_yield_strength'] * 1E-6,
        arguments['pressurant_tank_yield_strength'] * 1E-6
    ]

    ep_iter_col1 = [
        ep_cycle.combustion_chamber_pressure * 1E-6,
        ep_cycle.main_oxidizer_flow,
        ep_cycle.main_fuel_flow,
        ep_cycle.pumps_power_required * 1E-3,
        ep_cycle.characteristic_velocity,
        ep_cycle.ideal_thrust_coefficient
    ]
    ep_iter_col2 = [
        ep_cycle.battery.total_energy,
        ep_cycle.battery.energy_heat_loss,
        0
    ]
    ep_cycle.iterate_coolant_flow()
    ep_iter_col1 += [
        ep_cycle.oxidizer_pump.mass_flow,
        ep_cycle.fuel_pump.mass_flow,
        ep_cycle.oxidizer_pump.power_required * 1E-3,
        ep_cycle.fuel_pump.power_required * 1E-3,
        ep_cycle.pumps_power_required * 1E-3,
        ep_cycle.actual_battery_coolant_flow,
        ep_cycle.chamber_mass_flow
    ]

    ep_iter_col2 += [
        ep_cycle.actual_battery_coolant_flow,
        ep_cycle.oxidizer_pump.mass_flow,
        ep_cycle.fuel_pump.mass_flow,
        ep_cycle.oxidizer_pump.power_required * 1E-3,
        ep_cycle.fuel_pump.power_required * 1E-3,
        ep_cycle.pumps_power_required * 1E-3,
        ep_cycle.battery.total_energy,
        ep_cycle.battery.energy_heat_loss,
        ep_cycle.battery_coolant_temperature_change,
        None
    ]

    gg_before_iter_list = [
        gg_cycle.gas_generator.outlet_mass_flow,
        gg_cycle.gas_generator.mass_mixture_ratio / (
                gg_cycle.gas_generator.mass_mixture_ratio + 1) * gg_cycle.gas_generator.outlet_mass_flow,
        1 / (gg_cycle.gas_generator.mass_mixture_ratio + 1) * gg_cycle.gas_generator.outlet_mass_flow,
    ]

    gg_iter_list = [
        gg_cycle.combustion_chamber_pressure * 1E-6,
        gg_cycle.base_mass_flow * gg_cycle.mass_mixture_ratio / (gg_cycle.mass_mixture_ratio + 1),
        gg_cycle.base_mass_flow / (gg_cycle.mass_mixture_ratio + 1),
        gg_cycle.pumps_power_required * 1E-3,
        gg_cycle.characteristic_velocity,
        gg_cycle.ideal_thrust_coefficient
    ]

    gg_cycle.iterate_mass_flow()

    gg_iter_list += [
        gg_cycle.oxidizer_pump.mass_flow,
        gg_cycle.fuel_pump.mass_flow,
        gg_cycle.oxidizer_pump.power_required * 1E-3,
        gg_cycle.fuel_pump.power_required * 1E-3,
        gg_cycle.pumps_power_required * 1E-3,
        gg_cycle.chamber_mass_flow,
        gg_cycle.gg_mass_flow,
    ]

    gg_before_iter_list += [
        gg_cycle.oxidizer_pump.mass_flow,
        gg_cycle.fuel_pump.mass_flow,
        gg_cycle.oxidizer_pump.power_required * 1E-3,
        gg_cycle.fuel_pump.power_required * 1E-3,
        gg_cycle.pumps_power_required * 1E-3,
        None,
        gg_cycle.gas_generator.pressure * 1E-6,
        gg_cycle.oxidizer_pump.pressure_change * 1E-6,
        gg_cycle.fuel_pump.pressure_change * 1E-6,
        gg_cycle.gas_generator.gas_density
    ]

    gg_after_iter_list = [
        gg_cycle.gg_mass_flow,
        gg_cycle.gas_generator.mass_mixture_ratio / (
                gg_cycle.gas_generator.mass_mixture_ratio + 1) * gg_cycle.gas_generator.outlet_mass_flow,
        1 / (gg_cycle.gas_generator.mass_mixture_ratio + 1) * gg_cycle.gas_generator.outlet_mass_flow,
        gg_cycle.chamber_oxidizer_flow,
        gg_cycle.chamber_fuel_flow,
        gg_cycle.chamber_mass_flow,
        None,
        None,
        None,
        None,
        None,
        None,
        None
    ]

    gg_cycle_col = [
        gg_cycle.pumps_mass,
        gg_cycle.gas_generator.mass,
        gg_cycle.oxidizer.mass,
        gg_cycle.fuel.mass,
        gg_cycle.oxidizer.volume,
        gg_cycle.fuel.volume,
        gg_cycle.pressurant.mass,
        gg_cycle.pressurant_tank.volume,
        gg_cycle.oxidizer_tank.volume,
        gg_cycle.fuel_tank.volume,
        gg_cycle.oxidizer_tank.radius,
        gg_cycle.fuel_tank.radius,
        0.,
        0.,
        0.,
        gg_cycle.oxidizer_tank.initial_head,
        0.,
        0.,
        0.,
        gg_cycle.fuel_tank.initial_head,
        gg_cycle.oxidizer_tank.total_upper_pressure,
        gg_cycle.oxidizer_tank.total_lower_pressure,
        gg_cycle.fuel_tank.total_upper_pressure,
        gg_cycle.fuel_tank.total_lower_pressure,
        gg_cycle.oxidizer_tank.mass,
        gg_cycle.fuel_tank.mass,
        gg_cycle.pressurant_tank.mass,
        gg_cycle.mass,
        None,
        None,
        None
    ]
    ep_cycle_col = [
        ep_cycle.oxidizer_pump.mass,
        ep_cycle.fuel_pump.mass,
        ep_cycle.electric_motor.mass,
        ep_cycle.inverter.mass,
        ep_cycle.battery.mass,
        ep_cycle.oxidizer.mass,
        ep_cycle.fuel.mass,
        ep_cycle.oxidizer.volume,
        ep_cycle.fuel.volume,
        ep_cycle.pressurant.mass,
        ep_cycle.pressurant_tank.volume,
        ep_cycle.oxidizer_tank.volume,
        ep_cycle.fuel_tank.volume,
        ep_cycle.oxidizer_tank.radius,
        ep_cycle.fuel_tank.radius,
        0.,
        0.,
        0.,
        ep_cycle.oxidizer_tank.initial_head,
        0.,
        0.,
        0.,
        ep_cycle.fuel_tank.initial_head,
        ep_cycle.oxidizer_tank.total_upper_pressure,
        ep_cycle.oxidizer_tank.total_lower_pressure,
        ep_cycle.fuel_tank.total_upper_pressure,
        ep_cycle.fuel_tank.total_lower_pressure,
        ep_cycle.oxidizer_tank.mass,
        ep_cycle.fuel_tank.mass,
        ep_cycle.pressurant_tank.mass,
        ep_cycle.mass
    ]

    gg_list = [
        gg_cycle.feed_system_mass,
        gg_cycle.pumps_mass,
        gg_cycle.gas_generator.mass,
        0.,
        0.,
        0.,
        gg_cycle.tanks_mass,
        gg_cycle.oxidizer_tank.mass,
        gg_cycle.fuel_tank.mass,
        gg_cycle.pressurant_tank.mass,
        gg_cycle.pressurant.mass,
        gg_cycle.cc_propellant_mass + gg_cycle.gg_propellant_mass,
        gg_cycle.cc_propellant_mass,
        gg_cycle.gg_propellant_mass,
        0.0,
        gg_cycle.mass
    ]
    ep_list = [
        ep_cycle.feed_system_mass,
        0.,
        0.,
        ep_cycle.pumps_mass,
        ep_cycle.electric_motor.mass,
        ep_cycle.inverter.mass,
        ep_cycle.tanks_mass,
        ep_cycle.oxidizer_tank.mass,
        ep_cycle.fuel_tank.mass,
        ep_cycle.pressurant_tank.mass,
        ep_cycle.pressurant.mass,
        ep_cycle.props_mass + ep_cycle.battery.mass,
        ep_cycle.props_mass,
        0.0,
        ep_cycle.battery.mass,
        ep_cycle.mass
    ]

    gg_results = [
        gg_cycle.mass,
        gg_cycle.final_mass,
        gg_cycle.props_mass,
        gg_cycle.gg_propellant_mass,
        0.0,
        0.,
        gg_cycle.feed_system_mass / gg_cycle.props_mass,
        gg_cycle.mass_ratio,
        0.,
        gg_cycle.overall_specific_impulse,
        0.,
        4550,
        0.,
        gg_cycle.get_payload(4550),
    ]

    ep_results = [
        ep_cycle.mass,
        ep_cycle.final_mass,
        ep_cycle.props_mass,
        0.,
        ep_cycle.battery.mass,
        0.,
        ep_cycle.feed_system_mass / ep_cycle.props_mass,
        ep_cycle.mass_ratio,
        0.,
        ep_cycle.overall_specific_impulse,
        0.,
        4550,
        0.,
        ep_cycle.get_payload(4550),
    ]
    ratio_results = []
    for gg, ep in zip(gg_results, ep_results):
        new = gg/ep if ep!= 0. else 0.
        ratio_results.append(new)

    vals = {}
    for name in ['characteristic_velocity', 'ideal_thrust_coefficient', 'base_mass_flow', 'exit_pressure']:
        epval = getattr(ep_cycle, name)
        ggval = getattr(gg_cycle, name)
        if epval == ggval:
            vals[name] = epval
        else:
            print(f'{name} ep:{epval}, gg:{ggval}')

    col1 = input_col1 + [None] * 3 + gg_cycle_col
    empty_col = [None] * 70
    col2 = input_col2
    col3 = [None] * 9 + [vals['characteristic_velocity']] + [None] * 7 + gg_iter_list + [None] * 2 + ep_iter_col1 + [None] * 4 + gg_results
    colx = [None] * 9 + [vals['ideal_thrust_coefficient']] + [None] * 39 + ep_results
    col4 = [None] * 9 + [vals['base_mass_flow']] + [None] * 7 + gg_before_iter_list + [None] * 2 + ep_iter_col2 + [None]*4 + ratio_results
    col5 = [None] * 17 + gg_after_iter_list + [None] * 19 + gg_list
    col6 = [None] * 49 + ep_list
    col7 = [None] * 26 + ep_cycle_col
    coly = [None] * 9 + [vals['exit_pressure'] * 1e-2] + [None] * 30
    complete_sheet = pd.DataFrame(
        data=zip_longest(col1, empty_col, empty_col, col2, empty_col, empty_col, col3, colx, col4, empty_col, col5,
                         col6, empty_col, col7, coly, fillvalue=None))

    columns = ['GG', 'EP']
    results = pd.DataFrame(data=zip(gg_list, ep_list), columns=columns)
    inputs = pd.DataFrame(data=zip(input_col1, input_col2))
    details = pd.DataFrame(data=zip(gg_cycle_col, ep_cycle_col))
    gg_iters = pd.DataFrame(
        data=zip(gg_iter_list, [None] * len(gg_iter_list), gg_before_iter_list, [None] * len(gg_iter_list),
                 gg_after_iter_list))
    ep_iters = pd.DataFrame(data=zip(ep_iter_col1, [None] * len(ep_iter_col1), ep_iter_col2))
    suffix = 'Fixed' if kwak_fix else 'Normal'
    time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    with pd.ExcelWriter(rf'data\full_output_kwak_comparison{time}_{suffix}.xlsx') as writer:  # doctest: +SKIP
        for i, data in enumerate([results, inputs, details, gg_iters, ep_iters, complete_sheet]):
            data.to_excel(writer, sheet_name=f'Sheet{i}')

    rows = [
        'Feed System', 'Turbopump', 'GG', 'Pumps', 'Electric Motor', 'Inverter', 'Tanks', 'Oxidizer Tank', 'Fuel Tank',
        'Pressurant Tank', 'Pressurant', 'Props+Battery', 'CC Propellants', 'GG Propellants', 'Battery Pack', 'Total'
    ]
    results.index = rows
    pd.options.display.float_format = '{:.2f}'.format
    print(f'isp ep:{ep_cycle.overall_specific_impulse}, gg:{gg_cycle.overall_specific_impulse}')
    return results


def get_outputs():
    full_outputs = {}
    base_arguments2 = args.base_arguments.copy()
    ep_arguments2 = args.ep_arguments.copy()
    for burn_time, batt_coolant_flow in zip(burn_times, batt_coolant_flows):
        base_arguments2['burn_time'] = burn_time
        df = full_output(kwak_fix=True, design_args=args.desgin_arguments,
                         base_args=base_arguments2, ep_args=ep_arguments2)
        full_outputs[burn_time] = df
    return full_outputs


def get_compares():
    full_outputs = get_outputs()
    data = pd.read_excel('data/Kwak_Initial_Mass_Budgets.xlsx', index_col=0)
    data = data.fillna(0.0)
    kwak_values = {t_b: data[[f'GG {t_b}', f'EP {t_b}']] for t_b in burn_times}
    compares = {}
    for t_b in burn_times:
        compare_df = pd.concat([kwak_values[t_b], full_outputs[t_b]], axis=1)
        for name in ('GG', 'EP'):
            compare_df[f'{name} Diff.'] = compare_df[name] - compare_df[f'{name} {t_b}']
            compare_df[f'{name} Diff. %'] = compare_df[f'{name} Diff.'] / compare_df[f'{name} {t_b}'] * 100
        columns = list(compare_df.columns.values)
        rearrange = list(columns[:-3])
        rearrange.extend((columns[-2], columns[-3], columns[-1]))
        compare_df = compare_df[rearrange].fillna(0.0)
        compares[t_b] = compare_df
    return compares


def get_average_diff():
    compares = get_compares()
    average_diff_df = pd.DataFrame()
    for name in ('GG', 'EP'):
        average_diff_df[f'Average {name} Diff. %'] = (compares[300][f'{name} Diff. %']
                                                      + compares[390][f'{name} Diff. %']
                                                      + compares[1200][f'{name} Diff. %']) / 3
    return average_diff_df


def compares_to_csv(version: str) -> None:
    compares = get_compares()
    for burn_time in compares:
        compares[burn_time].to_csv(f'data/{burn_time}_full_output_comparison_table{version}.csv')


if __name__ == '__main__':
    # compares = get_compares()
    # compares_to_csv('2')
    design_args = {
        'thrust': 75e3,
        'combustion_chamber_pressure': 7e6,
        'burn_time': 500,
        'is_frozen': False
    }
    for kwak_fix in [True, False]:
        print(full_output(kwak_fix=kwak_fix, design_args=design_args))
