from math import radians
from EngineComponents.Abstract.Material import KwakGasGeneratorMaterial, KwakPressurantTankMaterial, \
    KwakPropellantTankMaterial, Inconel600, NarloyZ
from EngineComponents.Abstract.FlowState import ManualFlowState
from EngineFunctions.BaseFunctions import copy_without

desgin_arguments = {
    'thrust': 75e3,
    'combustion_chamber_pressure': 7e6,
    'burn_time': 500,
}

base_arguments_kwak = {
    'oxidizer_name': 'LO2_NASA',
    'fuel_name': 'RP1_NASA',
    'max_acceleration': 4.5 * 9.80665,
    'mass_mixture_ratio': 2.45,
    'oxidizer_initial_pressure': .4E6,
    'fuel_initial_pressure': .25E6,
    'oxidizer_pump_efficiency': .66,
    'fuel_pump_efficiency': .61,
    'pressurant_name': 'Helium',
    'pressurant_initial_pressure': 27E6,
    'pressurant_final_pressure': 5E6,
    'pressurant_initial_temperature': 100,
    'pressurant_margin_factor': 1.1,
    'pressurant_tank_safety_factor': 1.2,
    'propellant_margin_factor': 1.01,
    'tanks_structural_factor': 2.5,
    'ullage_volume_factor': 1.08,
    'fuel_tank_material': KwakPropellantTankMaterial,
    'oxidizer_tank_material': KwakPropellantTankMaterial,
    'pressurant_tank_material': KwakPressurantTankMaterial,
}

base_arguments_own = {
    'is_frozen': True,
    'oxidizer_initial_temperature': 90.19,
    'combustion_chamber_material': NarloyZ,
    'injector_material': NarloyZ,
    'nozzle_material': Inconel600,
    'combustion_chamber_safety_factor': 1.5,
    'injector_safety_factor': 1.5,
    'nozzle_safety_factor': 1.5,
    'injector_pressure_drop_factor': .15,
    'convergent_half_angle': radians(30),
    'convergent_throat_bend_ratio': 0.8,
    'convergent_chamber_bend_ratio': 1.0,
    'divergent_throat_half_angle': radians(15),
    'maximum_wall_temperature': 850,
    'thrust_chamber_wall_emissivity': .8,
    'hot_gas_emissivity': .1,
    'cooling_pressure_drop_factor': .4,
    'specific_impulse_quality_factor': None,
    'shaft_mechanical_efficiency': 0.95,
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
    'exhaust_material': Inconel600,
    'exhaust_safety_factor': 1.5,
}
gg_arguments = open_arguments | {
    'gg_stay_time': 10E-3,
    'gg_structural_factor': 2.5,
    'gg_material': KwakGasGeneratorMaterial,
}

cb_arguments = open_arguments | {}

oe_arguments = open_arguments | {}

oe1_arguments = oe_arguments | {
    '_secondary_fuel_pump_pressure_factor_first_guess': .4, 'secondary_fuel_pump_efficiency': None,
}

ep_arguments = {
    # All values taken from H.D. Kwak 2018 - "Performance assessment of electrically driven pump-fed LOX/kerosene
    # cycle rocket engine: Comparison with gas generator cycle"
    'fuel_pump_specific_power': 15E3,
    'oxidizer_pump_specific_power': 20E3,
    'electric_motor_specific_power': 5.3E3,
    'inverter_specific_power': 60E3,
    'battery_specific_power': 6.95E3,
    'battery_specific_energy': 198 * 3600,
    'electric_motor_efficiency': .95,
    'inverter_efficiency': .85,
    'battery_structural_factor': 1.2,
    'battery_coolant_temperature_change': 40,
    'electric_motor_heat_loss_factor': 0.015,
    'electric_motor_magnet_temp_limit': 400,
    'electric_motor_ox_leak_factor': 0.005,
}

# Kwak Arguments
gg_arguments_rp1_kwak = gg_arguments | {
    'gg_base_flow_state': ManualFlowState(propellant_name='ExhaustGas',
                                          temperature=open_arguments['turbine_maximum_temperature'],
                                          pressure=None,
                                          mass_flow=None,
                                          type='combusted',
                                          _molar_mass=0.03033368339292229,
                                          _specific_heat_capacity=2024.7,
                                          _heat_capacity_ratio=1.16,
                                          _density=None),

    'gg_mass_mixture_ratio': 0.320,
}

ep_arguments_rp1_kwak = ep_arguments | {'battery_coolant_specific_heat_capacity': 2009, }

common_arguments_kwak = base_arguments | {
    '_ignore_cooling': True,
    'fuel_initial_temperature': 263.6,
    'specific_impulse_quality_factor': 1.0,
    'shaft_mechanical_efficiency': 1.0,
    'exit_pressure_forced': 0.002E6,
    # To get the as close as possible to density given by Kwak with his initial pressure of 2.5 bar
    'oxidizer_initial_temperature': 93.340,  # To get the same density as Kwak with his initial pressure of 4 bar
    '_fuel_pump_pressure_factor_first_guess': 1.55,
    '_oxidizer_pump_pressure_factor_first_guess': 1.15,
    'is_frozen': False,

}

kwak_specific_arguments = {
    'pressurant_heat_capacity_ratio': 1.667,
    'pressurant_molar_mass': 0.00399733779,
    'manual_oxidizer_density': 1126.1,
    'manual_fuel_density': 804.2,
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
    'maximum_wall_temperature': 850,
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
