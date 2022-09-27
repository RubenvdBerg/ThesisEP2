from scipy.constants import g
from math import radians


def copy_without(origin_dict, iterable_keys):
    copy_dict = origin_dict.copy()
    for key in iterable_keys:
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
    # 'oxidizer_density': 1126.1,
    # 'fuel_density': 804.2,
    'tanks_material_density': 2850,
    'pressurant_tank_material_density': 4430,
    'tanks_yield_strength': 250E6,
    'pressurant_tank_yield_strength': 1100E6,
}
base_arguments_own = {
    'fuel_initial_temperature': 263.6,  # To get the as close as possible to density given by Kwak with his initial pressure of 2.5 bar
    'oxidizer_initial_temperature': 93.340,  # To get the same density as Kwak with his initial pressure of 4 bar
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
    'divergent_throat_half_angle': radians(35),
    'divergent_exit_half_angle': radians(5),
    'nozzle_type': 'conical',
    'maximum_wall_temperature': 850,
    'thrust_chamber_wall_emissivity': .8,
    'hot_gas_emissivity': .1,
    'convective_coefficient_mode': 'Modified Bartz',
    # 'coolant_liquid_heat_capacity': 1,
    # 'coolant_gas_heat_capacity': 1,
    # 'coolant_heat_of_vaporization': 1,
    # 'coolant_molar_mass': 1,
    # 'coolant_boiling_temp_1_bar': 1,
}

base_arguments = base_arguments_kwak | base_arguments_own

base_arguments_o = copy_without(base_arguments, ['mass_mixture_ratio'])


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


open_arguments = {
    'turbine_gas_specific_heat_capacity': 2024.7, 'turbine_gas_heat_capacity_ratio': 1.16,
    'turbine_pressure_ratio': 27, 'turbopump_specific_power': 13.5E3, 'turbine_efficiency': .52,
    'exhaust_thrust_contribution': .01,  'exhaust_expansion_ratio': 20
}
gg_arguments = open_arguments | {
    'turbine_maximum_temperature': 900, 'gg_mass_mixture_ratio': 0.320, 'gg_gas_specific_gas_constant': 274.1,
    'gg_stay_time': 10E-3, 'gg_structural_factor': 2.5,
    'gg_material_density': 8220, 'gg_yield_strength': 550E6
}

oe_arguments = open_arguments | {

}

ep_arguments = {
    'fuel_pump_specific_power': 15E3, 'oxidizer_pump_specific_power': 20E3, 'fuel_specific_heat': 2009,
    'electric_motor_specific_power': 5.3E3, 'inverter_specific_power': 60E3, 'battery_specific_power': 6.95E3,
    'battery_specific_energy': 198 * 3600, 'electric_motor_efficiency': .95, 'inverter_efficiency': .85,
    'battery_structural_factor': 1.2, 'battery_coolant_temperature_change': 40,
}

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
    'is_frozen': True,
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
    'is_frozen': True,
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
    'is_frozen': True,
    'exit_pressure_forced': None,
    'expansion_ratio_end_cooling': 15
}

se_21d_kwargs = base_arguments_o | {
    'thrust': 1947e3,
    'combustion_chamber_pressure': 6.649e6,
    'fuel_initial_temperature': 21,
    'oxidizer_initial_temperature': 90,
    'expansion_ratio': 12.52,
    'mass_mixture_ratio': 5.5,
    'fuel_name': 'LH2_NASA',
    # 'fuel_density': 70.22,
    'burn_time': 100,
    'is_frozen': True,
    'exit_pressure_forced': None,
    'expansion_ratio_end_cooling': 5,
    'chamber_characteristic_length': 4.0,
    'fuel_pump_efficiency': .7,
    'oxidizer_pump_efficiency': .76,
    'fuel_initial_pressure': .3e6,
    'oxidizer_initial_pressure': .5e6,
    'turbine_pressure_ratio': 27.7033333,
    'turbine_gas_specific_heat_capacity': None,
    'turbine_gas_heat_capacity_ratio': None,
    'turbopump_specific_power': 13.5E3,
    'turbine_efficiency': .45,
    'exhaust_thrust_contribution': .0084,
    'area_ratio_chamber_throat': (.985/2)**2/.286**2,
    'exhaust_expansion_ratio': 1.655,
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
    'is_frozen': True,
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
    'exhaust_thrust_contribution': .012,
    'area_ratio_chamber_throat': (.985/2)**2/.286**2,
    'divergent_throat_half_angle': radians(79),
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
    'is_frozen': True,
    'exit_pressure_forced': None,
    'expansion_ratio_end_cooling': None,
    'nozzle_type': 'conical',
    'maximum_wall_temperature': 850,
}

hyprob_kwargs = base_arguments_o | {
    'burn_time': 100,
    'throat_area': 0.002879753,
    'expansion_ratio': 8.848742189,
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
    'is_frozen': True,
    'thrust': 30e3,
}

le5a_kwargs_cnozzle = change_to_conical_nozzle(le5a_kwargs)

duel_pump_kwargs = {'fuel_pump_specific_power': 15E3, 'oxidizer_pump_specific_power': 20E3}

single_pump_kwargs = {'turbopump_specific_power': 13.5E3}
