import EngineArguments.arguments as args
from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle_DoublePump, OpenExpanderCycle
from KwakFix.KwakFixCycles import KwakEngineCycle


def get_default_kwargs(EngineClass: EngineCycle, mass_mixture_ratio: bool = True):
    if issubclass(EngineClass, KwakEngineCycle):
        default_args = args.common_arguments_kwak | args.kwak_specific_arguments
        if issubclass(EngineClass, ElectricPumpCycle):
            return default_args | args.ep_arguments_rp1_kwak
        else:
            return default_args | args.gg_arguments_rp1_kwak
    else:
        default_args = args.base_arguments if mass_mixture_ratio else args.base_arguments_o
        if issubclass(EngineClass, ElectricPumpCycle):
            return default_args | args.ep_arguments
        elif issubclass(EngineClass, GasGeneratorCycle):
            return default_args | args.gg_arguments
        elif issubclass(EngineClass, OpenExpanderCycle):
            return default_args | args.oe_arguments
        elif issubclass(EngineClass, OpenExpanderCycle_DoublePump):
            return default_args | args.oe1_arguments

if __name__ == '__main__':
    ep_kwargs = get_default_kwargs(ElectricPumpCycle)
    gg_kwargs = get_default_kwargs(GasGeneratorCycle)
    for key, val in get_default_kwargs(OpenExpanderCycle).items():
        if key not in ep_kwargs and key not in gg_kwargs:
            print(f"{key.replace('_',' ')},",val)