base_arguments = {
    'oxidizer_name': 'LO2_NASA', 'fuel_name': 'RP1_NASA',
    'exit_pressure': .002E6, 'max_acceleration': 4.5, 'heat_ratio_pressurant': 1.667,
    'mass_mixture_ratio': 2.45, 'pressurant_initial_pressure': 27E6, 'pressurant_final_pressure': 5E6,
    'oxidizer_initial_pressure': .4E6, 'fuel_initial_pressure': .25E6, 'fuel_pump_pressure_factor': 1.55,
    'oxidizer_pump_pressure_factor': 1.15, 'pressurant_gas_constant': 2080,
    'pressurant_initial_temperature': 100, 'oxidizer_pump_efficiency': .66, 'fuel_pump_efficiency': .61,
    'pressurant_margin_factor': 1.1, 'pressurant_tank_structural_factor': 1.2,
    'propellant_margin_factor': 1.01, 'tanks_structural_factor': 2.5, 'ullage_volume_factor': 1.08,
    'oxidizer_density': 1126.1, 'fuel_density': 804.2, 'tanks_material_density': 2850,
    'pressurant_tank_material_density': 4430, 'tanks_yield_strength': 250E6,
    'pressurant_tank_yield_strength': 1100E6
}
gg_arguments = {
    'gg_gas_specific_heat': 2024.7, 'heat_ratio_gg_gas': 1.16, 'mass_mixture_ratio_gg': 0.320,
    'turbine_pressure_ratio': 27, 'gas_constant_gg_gas': 274.1, 'turbine_inlet_temperature': 900,
    'gas_generator_stay_time': 10E-3, 'turbopump_specific_power': 13.5E3, 'turbine_efficiency': .52,
    'gg_structural_factor': 2.5, 'gg_material_density': 8220, 'gg_yield_strength': 550E6, 'gg_thrust_contribution': .01
}
ep_arguments = {
    'fuel_pump_specific_power': 15E3, 'oxidizer_pump_specific_power': 20E3, 'fuel_specific_heat': 2009,
    'electric_motor_specific_power': 5.3E3, 'inverter_specific_power': 60E3, 'battery_specific_power': 6.95E3,
    'battery_specific_energy': 198 * 3600, 'electric_motor_efficiency': .95, 'inverter_efficiency': .85,
    'battery_structural_factor': 1.2, 'coolant_allowable_temperature_change': 40,
}