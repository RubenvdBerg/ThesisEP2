from EngineComponents.Abstract.Material import StainlessSteel301_Annealed, Al6061T6
from EngineArguments.arguments import base_arguments, ep_arguments
base = base_arguments | ep_arguments
remove_keys = ['battery_coolant_temperature_change','electric_motor_heat_loss_factor', 'electric_motor_magnet_temp_limit', 'electric_motor_ox_leak_factor']
simple_base = {key:value for key,value in base.items() if key not in remove_keys}

lee_kwargs = simple_base | {
    'thrust': 500,
    'combustion_chamber_pressure': 20e5,
    'burn_time': 600,
    'ambient_pressure': 1028,
    'exit_pressure_forced': 1028,
    'tanks_structural_factor': 1,
    'propellant_margin_factor': 1.1,
    'injector_safety_factor': 2,
    'nozzle_safety_factor':2,
    'combustion_chamber_safety_factor':2,
    'max_acceleration': 0,
    'specific_impulse_quality_factor': 1,
    'combustion_chamber_material': StainlessSteel301_Annealed,
    'nozzle_material': StainlessSteel301_Annealed,
    'chamber_characteristic_length': 1.145,
    'fuel_tank_material': Al6061T6,
    'oxidizer_tank_material': Al6061T6,
    'oxidizer_initial_pressure': 3.79e5,
    'fuel_initial_pressure': 3.1e5,
    'ullage_volume_factor': 1.2,
    'fuel_pump_efficiency': .61,
    'oxidizer_pump_efficiency': .66,
    'electric_motor_efficiency': .87,
    'inverter_efficiency': .85,
    'electric_motor_specific_power': 875,
    'inverter_specific_power': 60e3,
    'battery_specific_power': 650,
    'battery_specific_energy': 325 * 3600,
    'battery_structural_factor': 1.2,
    '_ignore_cooling': True,
}

rutherford_kwargs = base | {
    'thrust': 24.91e3,
    'combustion_chamber_pressure': 7e6,
    'burn_time': 150,
    'expansion_ratio': 15,
    'mass_mixture_ratio': 2.4,
    'ambient_pressure': 101325,
}