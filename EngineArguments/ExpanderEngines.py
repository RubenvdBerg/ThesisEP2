from math import radians

from EngineArguments import arguments as args
from EngineArguments.arguments import base_arguments_o
from EngineComponents.Abstract.Material import NarloyZ

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
le5_kwargs = {'fuel_name': 'LH2_NASA',
              'thrust': 103e3,
              'combustion_chamber_pressure': 3.65e6,
              'expansion_ratio': 140,
              'mass_mixture_ratio': 5.5,
              'area_ratio_chamber_throat': 3.11,
              'chamber_characteristic_length': 0.84,
              'turbine_maximum_temperature': 837,
              'gg_pressure': 2.63e6,
              'fuel_pump_outlet_pressure': 6.04,
              'oxidizer_pump_outlet_pressure': 5.19,
              'fuel_pump_efficiency': 61.2e-2,
              'oxidizer_pump_efficiency': 65.3e-2,
              'burn_time': 0,
              'turbine_pressure_ratio': 0,
              'turbine_efficiency': 0,
              }

# Data from Kakuma 2000, Aoki 2001
le5b_kwargs = args.base_arguments_o | {
    'thrust': 137e3,
    'fuel_name': 'LH2_NASA',
    'mass_mixture_ratio': 5.0,
    'expansion_ratio': 110,
    'exit_pressure_forced': None,
    'ambient_pressure': 0,
    'burn_time': 534,
    'expansion_ratio_end_cooling': 6,  # Visual approximation
    'area_ratio_chamber_throat': 3.11,  # From McHugh1995
    'combustion_chamber_pressure': 3.62e6,
    'turbine_maximum_temperature': 409,
    'oxidizer_turbine_pressure_ratio': 2.172,
    'fuel_turbine_pressure_ratio': 5.125,
    'fuel_initial_temperature': 21,
    'oxidizer_initial_temperature': 90,
    # Assumed (from SE21D)
    'divergent_throat_half_angle': radians(15),
    'specific_impulse_quality_factor': .99,
    'oxidizer_turbine_efficiency': .40,
    'fuel_turbine_efficiency': .40,
    'shaft_mechanical_efficiency': .95,
    'oxidizer_exhaust_exit_pressure_forced': .04e6,
    'oxidizer_secondary_specific_impulse_quality_factor': .98,
    'exhaust_material': NarloyZ,
    'fuel_pump_efficiency': .75,
    'oxidizer_pump_efficiency': .8,
    'fuel_initial_pressure': 4e5,
    'oxidizer_initial_pressure': 2e5,
    # 'fuel_initial_pressure': 3e5,
    # 'oxidizer_initial_pressure': 6e5,
    'turbopump_specific_power': 13.5E3,
    # 'fuel_pump_outlet_pressure_forced': 6.86e6,
    # 'oxidizer_pump_outlet_pressure_forced': 5.3e6,
}

se_21d_kwargs = args.base_arguments_o | {
    'thrust': 1947.03e3,
    'combustion_chamber_pressure': 6.649e6,
    'fuel_initial_temperature': 21,
    'oxidizer_initial_temperature': 90,
    'expansion_ratio': 12.52,
    'mass_mixture_ratio': 5.5,
    'fuel_name': 'LH2_NASA',
    'burn_time': 100,
    'expansion_ratio_end_cooling': 5,
    'chamber_characteristic_length': 4.0,
    'fuel_pump_efficiency': .7,
    'secondary_fuel_pump_efficiency': .75,
    'oxidizer_pump_efficiency': .76,
    'fuel_initial_pressure': .3e6,
    'oxidizer_initial_pressure': .5e6,
    'turbine_maximum_temperature': 506.452,
    'turbopump_specific_power': 13.5E3,
    'area_ratio_chamber_throat': (.985 / 2) ** 2 / .286 ** 2,
    'ambient_pressure': 101325,
    'divergent_throat_half_angle': radians(15),
    'specific_impulse_quality_factor': .99,
    'oxidizer_turbine_efficiency': .45,
    'fuel_turbine_efficiency': .45,
    'shaft_mechanical_efficiency': .99,
    'oxidizer_exhaust_exit_pressure_forced': .04e6,
    'fuel_exhaust_exit_pressure_forced': .04e6,
    'oxidizer_turbine_outlet_pressure_forced': .3e6,
    'fuel_turbine_outlet_pressure_forced': .3e6,
    'oxidizer_secondary_specific_impulse_quality_factor': .98,
    'fuel_secondary_specific_impulse_quality_factor': .98,
    'exhaust_material': NarloyZ,
}
se_21d_vac_kwargs = se_21d_kwargs | {'thrust': 2196.682e3, 'ambient_pressure': 0}

if __name__ == '__main__':
    from EngineArguments.default_arguments import get_default_kwargs
    from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
    base_kwargs = get_default_kwargs(OpenExpanderCycle)
    kwargss = [le5b_kwargs]
    names = ['le5b']
    for name, kwargs in zip(names, kwargss):
        print('\n',name)
        for key, val in kwargs.items():
            if key in base_kwargs:
                if val != base_kwargs[key]:
                    print(key, val)
            else:
                print(key, val)
