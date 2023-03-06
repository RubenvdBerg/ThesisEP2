import EngineArguments.arguments as args
from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle_DoublePump
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
        elif issubclass(EngineClass, OpenExpanderCycle_DoublePump):
            return default_args | args.oe_arguments
