import arguments as args
from EngineCycles.GasGeneratorCycle.GGCycle import GasGeneratorCycle

engine = GasGeneratorCycle(**args.lrb_kwargs, **args.gg_arguments)
print(engine.cooling_channel_section.outlet_temperature)