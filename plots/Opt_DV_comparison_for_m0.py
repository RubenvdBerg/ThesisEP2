from plots.Easy_Data import get_individual_comparison_excel, get_engine
from plots.KwakPlots.Results_Comparison_RP1 import engine_kwargs
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle

kwargs_switcher = {
    ElectricPumpCycle: {
        'battery_specific_energy': 100000000.0,
        'battery_specific_power': 1E20,
        'combustion_chamber_pressure': 4.8752844098031405E6,
        'mass_mixture_ratio': 2.3451441054818876,
    },
    GasGeneratorCycle: {
        'combustion_chamber_pressure': 6.129431881640637E6,
        'mass_mixture_ratio': 2.363302462260184,
    },
    OpenExpanderCycle: {
        'combustion_chamber_pressure': 5.416918169462596E6,
        'mass_mixture_ratio': 2.354290510530121,
    },
}

engine_kwargs['burn_time'] = 500
engines = []
for engine_cycle in [ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle]:
    local_kwargs = kwargs_switcher[engine_cycle]
    total_kwargs = engine_kwargs | local_kwargs
    engine = get_engine(engine_cycle, total_kwargs)
    engines.append(engine)

get_individual_comparison_excel(engines, 'Opt_DV_comparison_for_m0.xlsx')