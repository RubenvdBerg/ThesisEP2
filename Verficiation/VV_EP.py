from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from KwakFix.KwakFixCycles import KwakFixElectricPumpCycle
from EngineArguments.default_arguments import get_default_kwargs
from EngineArguments.ElectricPumpEngines import rutherford_kwargs, lee_kwargs
from plots.Imaging.performance_image import make_performance_schematic
from plots.Imaging.mass_image import make_mass_schematic

design_kwargs = {
        'thrust': 100e3,
        'combustion_chamber_pressure': 100e5,
        'burn_time': 300,
        'exit_pressure_forced': 0.002e6,
    }

for cycle in [KwakFixElectricPumpCycle, ElectricPumpCycle]:
    kwargs = get_default_kwargs(cycle) | design_kwargs | {'oxidizer_initial_temperature': 93.34}
    engine = cycle(**kwargs)
    make_performance_schematic(engine)
