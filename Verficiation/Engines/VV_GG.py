import arguments as args

hm7b_kwargs = {'fuel_name': 'LH2_NASA',
               'thrust': 62.2e3,
               'combustion_chamber_pressure': 3.6e6,
               'expansion_ratio': 82.9,
               'mass_mixture_ratio': 4.565,
               'area_ratio_chamber_throat': 2.78,
               'chamber_characteristic_length': 0.68,
               'turbine_maximum_temperature': 860,
               'gg_pressure': 2.3e6,
               'fuel_pump_outlet_pressure': 5.55,
               'oxidizer_pump_outlet_pressure': 5.02,
               'fuel_pump_efficiency': 60e-2,
               'oxidizer_pump_efficiency': 73e-2,
               'turbine_pressure_ratio': 16.7,
               'turbine_efficiency': .45,
               'burn_time': 731,
               }

hm60_kwargs = {'fuel_name': 'LH2_NASA',
               'thrust': 1025e3,
               'combustion_chamber_pressure': 10e6,
               'expansion_ratio': 45,
               'mass_mixture_ratio': 5.1,
               'area_ratio_chamber_throat': 2.99,
               'chamber_characteristic_length': 0.84,
               'turbine_maximum_temperature': 871,
               'gg_pressure': 8.5e6,
               'fuel_pump_outlet_pressure': 15.8,
               'oxidizer_pump_outlet_pressure': 13,
               'fuel_pump_efficiency': 73e-2,
               'oxidizer_pump_efficiency': 76e-2,
               'fuel_turbine_pressure_ratio': 17,
               'oxidizer_turbine_pressure_ratio': 13.6,
               'fuel_turbine_efficiency': .59,
               'oxidizer_turbine_efficiency': .27,
               'burn_time': 605,
               }

j2_kwargs = {'fuel_name': 'LH2_NASA',
             'thrust': 1023e3,
             'combustion_chamber_pressure': 5.4e6,
             'expansion_ratio': 27.5,
             'mass_mixture_ratio': 5.5,
             'area_ratio_chamber_throat': 1.58,
             'chamber_characteristic_length': 0.62,
             'turbine_maximum_temperature': 922,
             'gg_pressure': 4.7e6,
             # 'fuel_pump_outlet_pressure': 8.62,
             # 'oxidizer_pump_outlet_pressure': 7.64,
             'fuel_pump_efficiency': 73e-2,
             'oxidizer_pump_efficiency': 80e-2,
             'fuel_turbine_pressure_ratio': 7.2,
             'oxidizer_turbine_pressure_ratio': 2.65,
             'fuel_turbine_efficiency': .6,
             'oxidizer_turbine_efficiency': .47,
             'burn_time': 475,
             }

# s4_kwarg = {'fuel_name': 'RP1_NASA',
#             'thrust': 364e3,
#             'combustion_chamber_pressure': 4.6e6,
#             'expansion_ratio': 25,
#             'mass_mixture_ratio': 2.27,
#             'area_ratio_chamber_throat': 1.66,
#             'chamber_characteristic_length': 1.09,
#             'turbine_maximum_temperature': 843.8,
#             'gg_pressure': 5.15e6,
#             'fuel_pump_outlet_pressure': 7.05,
#             'oxidizer_pump_outlet_pressure': 6.8,
#             'fuel_pump_efficiency': None,
#             'oxidizer_pump_efficiency': None,
#             'burn_time': None,
#             'turbine_pressure_ratio': None,
#             'turbine_efficiency': None,
#             }

h1_kwargs = {'fuel_name': 'RP1_NASA',
             'thrust': 945.4e3,
             'combustion_chamber_pressure': 4.12e6,
             'expansion_ratio': 8,
             'mass_mixture_ratio': 2.26,
             'area_ratio_chamber_throat': 1.67,
             'chamber_characteristic_length': 0.983,
             'turbine_maximum_temperature': 922,
             'gg_pressure': 4.22e6,
             'fuel_pump_outlet_pressure': 7.1,
             'oxidizer_pump_outlet_pressure': 6.3,
             'fuel_pump_efficiency': 71e-2,
             'oxidizer_pump_efficiency': 75e-2,
             'turbine_pressure_ratio': 18.21,
             'turbine_efficiency': .66,
             'burn_time': 150,
             }

rs27_kwargs = {'fuel_name': 'RP1_NASA',
               'thrust': 1043e3,
               'combustion_chamber_pressure': 4.87e6,
               'expansion_ratio': 12,
               'mass_mixture_ratio': 2.245,
               'area_ratio_chamber_throat': 1.62,
               'chamber_characteristic_length': 0.99,
               'turbine_maximum_temperature': 916,
               'gg_pressure': 4.7e6,
               'fuel_pump_outlet_pressure': 7.09,
               'oxidizer_pump_outlet_pressure': 7.25,
               'fuel_pump_efficiency': 71.8e-2,
               'oxidizer_pump_efficiency': 77.9e-2,
               'turbine_pressure_ratio': 221,
               'turbine_efficiency': .589,
               'burn_time': 274,
               }

f1_kwargs = {'fuel_name': 'RP1_NASA',
             'thrust': 7775.5e3,
             'combustion_chamber_pressure': 7.76e6,
             'expansion_ratio': 16,
             'mass_mixture_ratio': 2.27,
             'area_ratio_chamber_throat': "-",
             'chamber_characteristic_length': 1.22,
             'turbine_maximum_temperature': 1062,
             'gg_pressure': 6.76e6,
             'fuel_pump_outlet_pressure': 13,
             'oxidizer_pump_outlet_pressure': 11,
             'fuel_pump_efficiency': .714,  # Estimated from average!!
             'oxidizer_pump_efficiency': .764,  # Estimated from average!!
             'burn_time': 161,
             'turbine_pressure_ratio': 16.3,
             'turbine_efficiency': 0.605,
             }

gg_struct_kwargs = {
    'turbopump_specific_power': 13.5E3,
    'gg_stay_time': 10e-3,
    'gg_structural_factor': 2.5,
    'gg_material_density': 8220,
    'gg_yield_strength': 550E6
}

base_kwargs = args.base_arguments_o | {
    'ambient_pressure': 0,
    'exhaust_expansion_ratio': 4,
    'fuel_pump_efficiency': .70,
    'oxidizer_pump_efficiency': .70,
}

if __name__ == '__main__':
    from EngineCycles.BaseOpenCycle.OpenCycle import OpenEngineCycle_DoubleTurbine
    from EngineCycles.GasGeneratorCycle.GGCycle import GasGeneratorCycle, GasGeneratorCycle_DoubleTurbine
    from plots.Imaging.performance_image import make_schematic


    j2_total_kwargs = args.base_arguments_o | j2_kwargs | {
        'ambient_pressure': 0,
        'fuel_exhaust_expansion_ratio': 4,
        'oxidizer_exhaust_expansion_ratio': 4,
        'convergent_throat_bend_ratio': 0.4,
        'convergent_chamber_bend_ratio': 0.5,
    }
    engine = GasGeneratorCycle_DoubleTurbine(**j2_total_kwargs)
    make_schematic(engine)