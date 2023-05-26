from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineArguments.default_arguments import get_default_kwargs
from plots.KwakPlots.Results_Comparison_RP1_Check import compare_engine_masses


def get_engine(EngineClass: EngineCycle, **kwargs):
    default_kwargs = get_default_kwargs(EngineClass)
    total_kwargs = default_kwargs | kwargs
    engine = EngineClass(**total_kwargs)
    # make_performance_schematic(engine)
    return engine


if __name__ == '__main__':
    kwargs = {
        'is_frozen': True,
        'burn_time': 1200,
        'thrust': 100e3,
        'combustion_chamber_pressure': 3e6,
        'fuel_name': 'RP1_NASA',
        'mass_mixture_ratio': 2.45,
        'exit_pressure_forced': .002e6,
        'expansion_ratio_end_cooling': 10,
        'verbose': True,
        'maximum_wall_temperature': 900,
    }
    ep = get_engine(ElectricPumpCycle, **kwargs)
    gg = get_engine(GasGeneratorCycle, **kwargs)
    oe = get_engine(OpenExpanderCycle, **kwargs)
    compare_engine_masses(ep, gg, oe)

