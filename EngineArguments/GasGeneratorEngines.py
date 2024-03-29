import EngineArguments.arguments as args
from EngineComponents.Abstract.Material import Inconel600

gg_base_kwargs = args.base_arguments_o | {
    'ambient_pressure': 0,
    'fuel_pump_efficiency': .714,  # Estimated from average and overriden if known
    'oxidizer_pump_efficiency': .764,  # Estimated from average!!
    'turbopump_specific_power': 13.5E3,
    'gg_stay_time': args.gg_arguments['gg_stay_time'],
    'gg_structural_factor': args.gg_arguments['gg_structural_factor'],
    'gg_material': args.gg_arguments['gg_material'],
    'exhaust_material': Inconel600,
    'exhaust_safety_factor': 1.5,
}
exhaust_expansion_ratio = 4  # Arbitrarily Determined

# All values below taken from:
# McHugh 1995 - "Numerical Analysis of Existing Liquid Rocket Engines as a Design Process Starter"

# Double Turbine
hm60_kwargs = gg_base_kwargs | {
    'fuel_name': 'LH2_NASA',
    'thrust': 1025e3,
    'combustion_chamber_pressure': 10e6,
    'expansion_ratio': 45,
    'mass_mixture_ratio': 5.1,
    'area_ratio_chamber_throat': 2.99,
    'chamber_characteristic_length': 0.84,
    'turbine_maximum_temperature': 871,
    'gg_pressure': 8.5e6,
    'gg_mass_mixture_ratio': .9,
    # 'fuel_pump_outlet_pressure': 15.8,
    # 'oxidizer_pump_outlet_pressure': 13,
    'fuel_pump_efficiency': 73e-2,
    'oxidizer_pump_efficiency': 76e-2,
    'fuel_turbine_pressure_ratio': 17,
    'oxidizer_turbine_pressure_ratio': 13.6,
    'fuel_turbine_efficiency': .59,
    'oxidizer_turbine_efficiency': .27,
    'burn_time': 605,
    'fuel_exhaust_expansion_ratio': exhaust_expansion_ratio,
    'oxidizer_exhaust_expansion_ratio': exhaust_expansion_ratio,
}

j2_kwargs = gg_base_kwargs | {
    'fuel_name': 'LH2_NASA',
    'thrust': 1023e3,
    'combustion_chamber_pressure': 5.4e6,
    'expansion_ratio': 27.5,
    'mass_mixture_ratio': 5.5,
    'area_ratio_chamber_throat': 1.58,
    'convergent_throat_bend_ratio': 0.4,
    'convergent_chamber_bend_ratio': 0.5,
    'chamber_characteristic_length': 0.62,
    'turbine_maximum_temperature': 922,
    'gg_pressure': 4.7e6,
    'gg_mass_mixture_ratio': .94,
    # 'fuel_pump_outlet_pressure': 8.62,
    # 'oxidizer_pump_outlet_pressure': 7.64,
    'fuel_pump_efficiency': 73e-2,
    'oxidizer_pump_efficiency': 80e-2,
    'fuel_turbine_pressure_ratio': 7.2,
    'oxidizer_turbine_pressure_ratio': 2.65,
    'fuel_turbine_efficiency': .6,
    'oxidizer_turbine_efficiency': .47,
    'burn_time': 475,
    'oxidizer_exhaust_expansion_ratio': exhaust_expansion_ratio,
}

# Single Turbine
hm7b_kwargs = gg_base_kwargs | {
    'fuel_name': 'LH2_NASA',
    'thrust': 62.2e3,
    'combustion_chamber_pressure': 3.6e6,
    'expansion_ratio': 82.9,
    'mass_mixture_ratio': 4.565,
    'area_ratio_chamber_throat': 2.78,
    'chamber_characteristic_length': 0.68,
    'turbine_maximum_temperature': 860,
    'gg_pressure': 2.3e6,
    'gg_mass_mixture_ratio': .87,
    # 'fuel_pump_outlet_pressure': 5.55,
    # 'oxidizer_pump_outlet_pressure': 5.02,
    'fuel_pump_efficiency': 60e-2,
    'oxidizer_pump_efficiency': 73e-2,
    'turbine_pressure_ratio': 16.7,
    'turbine_efficiency': .45,
    'burn_time': 731,
    'exhaust_expansion_ratio': exhaust_expansion_ratio,
}

h1_kwargs = gg_base_kwargs | {
    'fuel_name': 'RP1_NASA',
    'thrust': 945.4e3,
    'combustion_chamber_pressure': 4.12e6,
    'expansion_ratio': 8,
    'mass_mixture_ratio': 2.26,
    'area_ratio_chamber_throat': 1.67,
    'chamber_characteristic_length': 0.983,
    'turbine_maximum_temperature': 922,
    'gg_pressure': 4.22e6,
    'gg_mass_mixture_ratio': .342,
    # 'fuel_pump_outlet_pressure': 7.1,
    # 'oxidizer_pump_outlet_pressure': 6.3,
    'fuel_pump_efficiency': 71e-2,
    'oxidizer_pump_efficiency': 75e-2,
    'turbine_pressure_ratio': 18.21,
    'turbine_efficiency': .66,
    'burn_time': 150,
    'exhaust_expansion_ratio': exhaust_expansion_ratio,
}

rs27_kwargs = gg_base_kwargs | {
    'fuel_name': 'RP1_NASA',
    'thrust': 1043e3,
    'combustion_chamber_pressure': 4.87e6,
    'expansion_ratio': 12,
    'mass_mixture_ratio': 2.245,
    'area_ratio_chamber_throat': 1.62,
    'convergent_throat_bend_ratio': 0.4,
    'convergent_chamber_bend_ratio': 0.5,
    'chamber_characteristic_length': 0.99,
    'turbine_maximum_temperature': 916,
    'gg_pressure': 4.7e6,
    'gg_mass_mixture_ratio': .33,
    # 'fuel_pump_outlet_pressure': 7.09,
    # 'oxidizer_pump_outlet_pressure': 7.25,
    'fuel_pump_efficiency': 71.8e-2,
    'oxidizer_pump_efficiency': 77.9e-2,
    'turbine_pressure_ratio': 22.0,
    'turbine_efficiency': .589,
    'burn_time': 274,
    'exhaust_expansion_ratio': exhaust_expansion_ratio,
}

f1_kwargs = gg_base_kwargs | {
    'fuel_name': 'RP1_NASA',
    'thrust': 7775.5e3,
    'combustion_chamber_pressure': 7.76e6,
    'expansion_ratio': 16,
    'mass_mixture_ratio': 2.27,
    'area_ratio_chamber_throat': None,
    'chamber_characteristic_length': 1.22,
    'turbine_maximum_temperature': 1062,
    'gg_pressure': 6.76e6,
    'gg_mass_mixture_ratio': .416,
    # 'fuel_pump_outlet_pressure': 13,
    # 'oxidizer_pump_outlet_pressure': 11,
    # 'fuel_pump_efficiency': None,
    # 'oxidizer_pump_efficiency': None,
    'burn_time': 161,
    'turbine_pressure_ratio': 16.3,
    'turbine_efficiency': 0.605,
    'exhaust_expansion_ratio': exhaust_expansion_ratio,
}

s4_kwarg = gg_base_kwargs | {
    'fuel_name': 'RP1_NASA',
    'thrust': 364e3,
    'combustion_chamber_pressure': 4.6e6,
    'expansion_ratio': 25,
    'mass_mixture_ratio': 2.27,
    'area_ratio_chamber_throat': 1.66,
    'chamber_characteristic_length': 1.09,
    'turbine_maximum_temperature': 843.8,
    'gg_pressure': 5.15e6,
    'fuel_pump_outlet_pressure': 7.05,
    'oxidizer_pump_outlet_pressure': 6.8,
    'fuel_pump_efficiency': None,
    'oxidizer_pump_efficiency': None,
    'burn_time': None,
    'turbine_pressure_ratio': None,
    'turbine_efficiency': None,
    'exhaust_expansion_ratio': exhaust_expansion_ratio,
}

