from scipy.constants import g
from math import radians


def copy_without(origin_dict, key):
    copy_dict = origin_dict.copy()
    copy_dict.pop(key)
    return copy_dict


desgin_arguments = {
    'thrust': 75e3,
    'combustion_chamber_pressure': 7e6,
    'burn_time': 500,
    'is_frozen': False
}

base_arguments_kwak = {
    'oxidizer_name': 'LO2_NASA',
    'fuel_name': 'RP1_NASA',
    'exit_pressure_forced': 0.002E6,
    'max_acceleration': 4.5 * 9.80665,
    'pressurant_heat_capacity_ratio': 1.667,
    'mass_mixture_ratio': 2.45,
    'pressurant_initial_pressure': 27E6,
    'pressurant_final_pressure': 5E6,
    'oxidizer_initial_pressure': .4E6,
    'fuel_initial_pressure': .25E6,
    'fuel_pump_pressure_factor': 1.55,
    'oxidizer_pump_pressure_factor': 1.15,
    'pressurant_molar_mass': 0.00399733779,  # From gas_constant 2048
    'pressurant_initial_temperature': 100,
    'oxidizer_pump_efficiency': .66,
    'fuel_pump_efficiency': .61,
    'pressurant_margin_factor': 1.1,
    'pressurant_tank_safety_factor': 1.2,
    'propellant_margin_factor': 1.01,
    'tanks_structural_factor': 2.5,
    'ullage_volume_factor': 1.08,
    'oxidizer_density': 1126.1,
    'fuel_density': 804.2,
    'tanks_material_density': 2850,
    'pressurant_tank_material_density': 4430,
    'tanks_yield_strength': 250E6,
    'pressurant_tank_yield_strength': 1100E6,
}
base_arguments_own = {
    'combustion_chamber_material_density': 8470,
    'combustion_chamber_yield_strength': 300e6,
    'combustion_chamber_safety_factor': 1,
    'injector_material_density': 8470,
    'injector_yield_strength': 300e6,
    'injector_safety_factor': 1,
    'injector_propellant_is_gas': False,
    'convergent_half_angle': radians(30),
    'convergent_throat_bend_ratio': 0.8,
    'convergent_chamber_bend_ratio': 1.0,
    'chamber_throat_area_ratio': (80 / 50) ** 2,
    'divergent_throat_half_angle': radians(35),
    'divergent_exit_half_angle': radians(5),
    'nozzle_type': 'bell',
    'maximum_wall_temperature': 850,
    'thrust_chamber_wall_emissivity': .8,
    'hot_gas_emissivity': .1,
    'convective_coefficient_mode': 'Modified Bartz',
    'coolant_liquid_heat_capacity': 1,
    'coolant_gas_heat_capacity': 1,
    'coolant_heat_of_vaporization': 1,
    'coolant_molar_mass': 1,
    'coolant_boiling_temp_1_bar': 1,
    'coolant_inlet_temperature': 1
}

base_arguments = base_arguments_kwak | base_arguments_own

base_arguments_o = copy_without(base_arguments, 'mass_mixture_ratio')

open_arguments = {
    'turbine_gas_specific_heat_capacity': 2024.7, 'turbine_gas_heat_capacity_ratio': 1.16,
    'turbine_pressure_ratio': 27, 'turbopump_specific_power': 13.5E3, 'turbine_efficiency': .52,

}
gg_arguments = open_arguments | {
    'turbine_maximum_temperature': 900, 'gg_mass_mixture_ratio': 0.320, 'gg_gas_gas_constant': 274.1, 'gg_stay_time': 10E-3,'gg_structural_factor': 2.5,
    'gg_material_density': 8220, 'gg_yield_strength': 550E6, 'exhaust_thrust_contribution': .01
}

oe_arguments = open_arguments | {
    'turbine_inlet_temperature': 900
}

ep_arguments = {
    'fuel_pump_specific_power': 15E3, 'oxidizer_pump_specific_power': 20E3, 'fuel_specific_heat': 2009,
    'electric_motor_specific_power': 5.3E3, 'inverter_specific_power': 60E3, 'battery_specific_power': 6.95E3,
    'battery_specific_energy': 198 * 3600, 'electric_motor_efficiency': .95, 'inverter_efficiency': .85,
    'battery_structural_factor': 1.2, 'coolant_allowable_temperature_change': 40,
}

liquid_oxygen_coolant = {
    'coolant_liquid_heat_capacity': 1,
    'coolant_gas_heat_capacity': 1,
    'coolant_heat_of_vaporization': float(213./15.999),
    'coolant_molar_mass': 15.999E-3,
    'coolant_boiling_temp_1_bar': 90.15,
    'coolant_inlet_temperature': 90.15
}

rp1_coolant = {
    'coolant_liquid_heat_capacity': 1,
    'coolant_gas_heat_capacity': None,
    'coolant_heat_of_vaporization': float(246./170.),
    'coolant_molar_mass': 170.0E-3,
    'coolant_boiling_temp_1_bar': 489.45,
    'coolant_inlet_temperature': 293.15
}

liquid_hydrogen_coolant = {
    'coolant_liquid_heat_capacity': 1,
    'coolant_gas_heat_capacity': 1,
    'coolant_heat_of_vaporization': 446./2.01588,
    'coolant_molar_mass': 2.01588E-3,
    'coolant_boiling_temp_1_bar': 20.25,
    'coolant_inlet_temperature': 20.25
}