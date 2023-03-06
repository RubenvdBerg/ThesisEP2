import EngineArguments.arguments as args
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from plots.Imaging.performance_image import make_performance_schematic

vac_ruth = {
    "thrust": 25800,
    "burn_time": 360,
    "combustion_chamber_pressure": 14e6,
    "expansion_ratio": 14,
}
sea_ruth = {
    "thrust": 24910,
    "burn_time": 150,
    "combustion_chamber_pressure": 14e6,
    "expansion_ratio": 14,
}

engine = ElectricPumpCycle(
    **sea_ruth,
    **args.base_arguments,
    **args.ep_arguments)
make_performance_schematic(engine)
