import json

from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineArguments.default_arguments import get_default_kwargs
from plots.KwakPlots.Results_Comparison_RP1 import engine_kwargs
from plots.Opt_DV_for_years import get_tiede_battery_specific_energy,get_sakama_electric_motor_specific_power

path = r'C:\Users\rvand\PycharmProjects\ThesisEP2\data\opt_compare\20230615-112801_DV_Opt_Comparison_YEARS_Ignore_DP'
with open(path, 'r') as file:
    dict = json.load(file)
info = dict.pop('info')
for tb in dict:
    burn_time = int(tb)
    local_dict = dict[tb]['ElectricPumpCycle']
    for i in range(2):
        years = info['input_vals'][i]
        if years == 2067.0:
            years = 2067.65
        def_kwargs = get_default_kwargs(ElectricPumpCycle, False)
        total_kwargs = def_kwargs | info['engine_inputs'] | {
            'battery_specific_energy': get_tiede_battery_specific_energy(years, 1),
            'electric_motor_specific_power': get_sakama_electric_motor_specific_power(years),
            'burn_time': burn_time,
            'battery_specific_power':1E20,
            'combustion_chamber_pressure': local_dict['chamber_pressure'][i]*1e6,
            'mass_mixture_ratio': local_dict['mixture_ratio'][i],
        }
        engine = ElectricPumpCycle(**total_kwargs, verbose=False)
        print(f'{burn_time:.0f} {years:.1f} {engine.battery.specific_power_required}')