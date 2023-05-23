from EngineArguments.default_arguments import get_default_kwargs
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from Results_Comparison_RP1 import engine_kwargs
from plots.Easy_Data import get_comparison_excel


def check_gg_vs_ep():
    engine_kwargs = {
        'thrust': 100e3,
        'burn_time': 300,
        'exit_pressure_forced': 0.002e6,
        'expansion_ratio_end_cooling': 10,
        'combustion_chamber_pressure': 10e6,
        'maximum_wall_temperature': 900,
    }
    get_comparison_excel('New_Test_Excel_GG_vs_EP',
                         input_attribute='combustion_chamber_pressure',
                         input_range=(3e6,10e6),
                         input_nums=3,
                         input_prefix='M',
                         **engine_kwargs)

if __name__ == '__main__':
    check_gg_vs_ep()