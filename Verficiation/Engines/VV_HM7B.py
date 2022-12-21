import arguments as args
from EngineCycles.GasGeneratorCycle.GGCycle import GasGeneratorCycle
from plots.Imaging.performance_image import make_schematic

hm7b_kwargs = args.base_arguments_o | {
    'burn_time': 950,
    'thrust': 62.6e3,
    'ambient_pressure': 0,
    'fuel_name': 'LH2_NASA',
    'expansion_ratio': 82.9,
    'fuel_initial_pressure': 2e5,
    'oxidizer_initial_pressure': 3e5,
    'fuel_pump_efficiency': .73,
    'oxidizer_pump_efficiency': .60,
    'combustion_chamber_pressure': 36e5,
    'mass_mixture_ratio': 5.15,
    'turbine_maximum_temperature': 860,
    'turbine_pressure_ratio': 16.7,
    'turbine_efficiency': .45,
    'exhaust_expansion_ratio': 4,

    # Required but unnecessary for flow estimate
    'turbopump_specific_power': 13.5E3,
    'gg_stay_time': 10e-3,
    'gg_structural_factor': 2.5,
    'gg_material_density': 8220,
    'gg_yield_strength': 550E6
}
hm7b_engine = GasGeneratorCycle(**hm7b_kwargs)
make_schematic(hm7b_engine)