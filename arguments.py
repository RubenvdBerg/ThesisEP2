from scipy.constants import g
from math import radians


def copy_without(origin_dict, iterable_keys):
    copy_dict = origin_dict.copy()
    for key in iterable_keys:
        copy_dict.pop(key)
    return copy_dict


def change_to_conical_nozzle(arg_dict, throat_half_angle=radians(15)):
    copy_dict = copy_without(arg_dict, [
        'divergent_throat_half_angle',
        'divergent_exit_half_angle',
        'nozzle_type'
    ])
    new_args_dict = {'divergent_throat_half_angle': throat_half_angle,
                     'divergent_exit_half_angle': None,
                     'nozzle_type': 'conical'}
    return copy_dict | new_args_dict


desgin_arguments = {
    'thrust': 75e3,
    'combustion_chamber_pressure': 7e6,
    'burn_time': 500,
    'is_frozen': False
}

base_arguments_kwak = {
    'oxidizer_name': 'LO2_NASA',
    'fuel_name': 'RP1_NASA',
    'max_acceleration': 4.5 * 9.80665,
    'pressurant_heat_capacity_ratio': 1.667,
    'mass_mixture_ratio': 2.45,
    'pressurant_initial_pressure': 27E6,
    'pressurant_final_pressure': 5E6,
    'oxidizer_initial_pressure': .4E6,
    'fuel_initial_pressure': .25E6,
    'pressurant_molar_mass': 0.00399733779,  # From gas_constant 2048
    'pressurant_initial_temperature': 100,
    'oxidizer_pump_efficiency': .66,
    'fuel_pump_efficiency': .61,
    'pressurant_margin_factor': 1.1,
    'pressurant_tank_safety_factor': 1.2,
    'propellant_margin_factor': 1.01,
    'tanks_structural_factor': 2.5,
    'ullage_volume_factor': 1.08,
    'tanks_material_density': 2850,
    'pressurant_tank_material_density': 4430,
    'tanks_yield_strength': 250E6,
    'pressurant_tank_yield_strength': 1100E6,
}

base_arguments_own = {
    'is_frozen': True,
    'oxidizer_initial_temperature': 90.19,
    'combustion_chamber_material_density': 8470,
    'combustion_chamber_yield_strength': 300e6,
    'combustion_chamber_safety_factor': 1,
    'injector_material_density': 8470,
    'injector_yield_strength': 300e6,
    'injector_safety_factor': 1,
    'injector_propellant_is_gas': False,
    'injector_pressure_drop_factor': .15,
    'convergent_half_angle': radians(30),
    'convergent_throat_bend_ratio': 0.8,
    'convergent_chamber_bend_ratio': 1.0,
    'divergent_throat_half_angle': radians(15),
    'nozzle_type': 'conical',
    'maximum_wall_temperature': 850,
    'thrust_chamber_wall_emissivity': .8,
    'hot_gas_emissivity': .1,
    'convective_coefficient_mode': 'Modified Bartz',
    'cooling_pressure_drop_factor': .4,
    'specific_impulse_correction_factor': 1.0,
    'shaft_mechanical_efficiency': 1.0,
}

duel_pump_kwargs = {'fuel_pump_specific_power': 15E3, 'oxidizer_pump_specific_power': 20E3}
single_pump_kwargs = {'turbopump_specific_power': 13.5E3}

base_arguments = base_arguments_kwak | base_arguments_own
base_arguments_o = copy_without(base_arguments, ['mass_mixture_ratio'])

open_arguments = {
    'turbine_pressure_ratio': 27,
    'turbine_efficiency': .52,
    'turbine_maximum_temperature': 900,
    'turbopump_specific_power': 13.5E3,
    'exhaust_expansion_ratio': 20,
}
gg_arguments = open_arguments | {
    'gg_stay_time': 10E-3,
    'gg_structural_factor': 2.5,
    'gg_material_density': 8220,
    'gg_yield_strength': 550E6
}

cb_arguments = open_arguments | {}

oe_arguments = open_arguments | {
    'secondary_fuel_pump_pressure_change_factor': .4, 'secondary_fuel_pump_efficiency': None,
}

ep_arguments = {
    'fuel_pump_specific_power': 15E3, 'oxidizer_pump_specific_power': 20E3,
    'electric_motor_specific_power': 5.3E3, 'inverter_specific_power': 60E3, 'battery_specific_power': 6.95E3,
    'battery_specific_energy': 198 * 3600, 'electric_motor_efficiency': .95, 'inverter_efficiency': .85,
    'battery_structural_factor': 1.2, 'battery_coolant_temperature_change': 40,
    'electric_motor_heat_loss_factor': 0.015,
    'electric_motor_magnet_temp_limit': 400,
    'electric_motor_ox_leak_factor': 0.005,
}

# Kwak Arguments
gg_arguments_rp1_kwak = gg_arguments | {
    'gg_gas_specific_heat_capacity': 2024.7,
    'gg_gas_heat_capacity_ratio': 1.16,
    'gg_gas_molar_mass': 0.03033368339292229,
    'gg_mass_mixture_ratio': 0.320,
}

ep_arguments_rp1_kwak = ep_arguments | {'battery_coolant_specific_heat_capacity': 2009, }

common_arguments_kwak = base_arguments | {
    '_ignore_cooling': True,
    'oxidizer_density': 1126.1,
    'fuel_density': 804.2,
    'fuel_initial_temperature': 263.6,
    # To get the as close as possible to density given by Kwak with his initial pressure of 2.5 bar
    'oxidizer_initial_temperature': 93.340,  # To get the same density as Kwak with his initial pressure of 4 bar
}

# Propellant properties dicts
liquid_oxygen_coolant = {
    'coolant_liquid_heat_capacity': 1,
    'coolant_gas_heat_capacity': 1,
    'coolant_heat_of_vaporization': float(213. / 15.999),
    'coolant_molar_mass': 15.999E-3,
    'coolant_boiling_temp_1_bar': 90.15,
    'coolant_inlet_temperature': 90.15
}

rp1_coolant = {
    'coolant_liquid_heat_capacity': 1,
    'coolant_gas_heat_capacity': None,
    'coolant_heat_of_vaporization': float(246. / 170.),
    'coolant_molar_mass': 170.0E-3,
    'coolant_boiling_temp_1_bar': 489.45,
    'coolant_inlet_temperature': 293.15
}

liquid_hydrogen_coolant = {
    'coolant_liquid_heat_capacity': 1,
    'coolant_gas_heat_capacity': 1,
    'coolant_heat_of_vaporization': 446. / 2.01588,
    'coolant_molar_mass': 2.01588E-3,
    'coolant_boiling_temp_1_bar': 20.25,
    'coolant_inlet_temperature': 20.25
}

# Engine arguments
tcd1_kwargs = base_arguments_o | {
    'fuel_initial_temperature': 20.25,
    'thrust': 100e3,
    'combustion_chamber_pressure': 55e5,
    'expansion_ratio': 143.2,
    'mass_mixture_ratio': 5.6,
    'area_ratio_chamber_throat': (80 / 55) ** 2,
    'chamber_characteristic_length': 1.095,
    'fuel_name': 'LH2_NASA',
    'burn_time': 100,
    'exit_pressure_forced': None,
    'expansion_ratio_end_cooling': 22,
    'nozzle_type': 'conical',
    'maximum_wall_temperature': 850,
}

le5a_kwargs = base_arguments_o | {
    'thrust': 121.3e3,
    'combustion_chamber_pressure': 40e5,
    'expansion_ratio': 130,
    'mass_mixture_ratio': 5,
    'fuel_name': 'LH2_NASA',
    'burn_time': 609,
    'exit_pressure_forced': None,
    'expansion_ratio_end_cooling': 30
}

lrb_kwargs = base_arguments_o | {
    'thrust': 2118.14e3,
    'combustion_chamber_pressure': 65e5,
    'expansion_ratio': 15,
    'mass_mixture_ratio': 2.4,
    'fuel_name': 'RP1_NASA',
    'burn_time': 100,
    'exit_pressure_forced': None,
    'expansion_ratio_end_cooling': 15
}

mira_kwargs = base_arguments_o | {
    'thrust': 100e3,
    'combustion_chamber_pressure': 48e5,
    'expansion_ratio': 125,
    'mass_mixture_ratio': 3.9,
    'fuel_name': 'CH4',
    # 'fuel_density': 423.9,
    'fuel_initial_temperature': 115,
    'oxidizer_initial_temperature': 90,
    'burn_time': 100,
    'exit_pressure_forced': None,
    'maximum_wall_temperature': 800,
    'distance_from_throat_start_cooling': None,
    'distance_from_throat_end_cooling': .4,
    'chamber_characteristic_length': 1.3,
    'fuel_pump_efficiency': .7,
    'oxidizer_pump_efficiency': .7,
    'fuel_initial_pressure': 3e5,
    'oxidizer_initial_pressure': 3e5,
    'fuel_pump_pressure_factor': 1.6,
    'oxidizer_pump_pressure_factor': 1.3,
    'turbine_pressure_ratio': 14.478260869565217391304347826087,
    'turbine_gas_specific_heat_capacity': None,
    'turbine_gas_heat_capacity_ratio': None,
    'turbopump_specific_power': 13.5E3,
    'turbine_efficiency': .6,
    'area_ratio_chamber_throat': (.985 / 2) ** 2 / .286 ** 2,
    'divergent_throat_half_angle': radians(25),
    'exhaust_expansion_ratio': 20

}

vinci_kwargs = base_arguments_o | {
    'fuel_initial_temperature': 20.25,
    'thrust': 180e3,
    'combustion_chamber_pressure': 60.8e5,
    'expansion_ratio': 240,
    'mass_mixture_ratio': 5.8,
    'area_ratio_chamber_throat': None,
    'chamber_characteristic_length': None,
    'fuel_name': 'LH2_NASA',
    'burn_time': 720,
    'exit_pressure_forced': None,
    'expansion_ratio_end_cooling': None,
    'nozzle_type': 'conical',
    'maximum_wall_temperature': 850,
}

hyprob_kwargs = base_arguments_o | {
    'burn_time': 100,
    'expansion_ratio': 8.848742189,
    'exit_pressure_forced': None,
    'area_ratio_chamber_throat': 3.849543342,
    'combustion_chamber_pressure': 56e5,
    'mass_mixture_ratio': 3.5,
    'fuel_initial_temperature': 110,
    'fuel_name': 'CH4',
    'distance_from_throat_start_cooling': None,
    'distance_from_throat_end_cooling': None,
    'chamber_characteristic_length': 0.901651725,
    'divergent_throat_half_angle': radians(19.88516511),
    'divergent_exit_half_angle': None,
    'convergent_half_angle': radians(23.81382356),
    'thrust': 30e3,
}

denies_kwargs = base_arguments_o | {
    'burn_time': 100,
    'expansion_ratio': 15.25037158,
    'exit_pressure_forced': None,
    'area_ratio_chamber_throat': 11.89060642,
    'combustion_chamber_pressure': 40e5,
    'mass_mixture_ratio': 3.16,
    'fuel_initial_temperature': 110,
    'fuel_name': 'CH4',
    'chamber_characteristic_length': 1.75,
    'divergent_throat_half_angle': radians(40),
    'divergent_exit_half_angle': None,
    'convergent_half_angle': radians(60),
    'thrust': 10e3,
    'convergent_throat_bend_ratio': 0.5,
    'convergent_chamber_bend_ratio': .4,
}

# le5a_kwargs_cnozzle = change_to_conical_nozzle(le5a_kwargs)
