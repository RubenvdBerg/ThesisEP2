from plots.single_engine_comparison import plot_mass_pie_chart_comparison
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
import EngineArguments.arguments as args
from plots.Imaging.mass_image import make_mass_schematic
from plots.Imaging.performance_image import make_performance_schematic

design_args = {'thrust': 1000e3,
               'burn_time': 900,
               'combustion_chamber_pressure': 3e6,
               'expansion_ratio': 60,
               'ambient_pressure': 0}

# plot_mass_pie_chart_comparison(design_args)

extra_args_switch = {
    # ElectricPumpCycle: args.ep_arguments,
    GasGeneratorCycle: args.gg_arguments,
    OpenExpanderCycle: args.oe_arguments,
}

for CycleClass, unique_args in extra_args_switch.items():
    total_args = args.base_arguments | unique_args | design_args
    engine = CycleClass(**total_args)
    print(engine.secondary_exhaust.inlet_flow_state.molar_mass)
    make_performance_schematic(engine)
    make_mass_schematic(engine)