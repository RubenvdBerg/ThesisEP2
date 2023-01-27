from KwakFix.KwakFixCycles import KwakFixElectricPumpCycle, KwakFixGasGeneratorCycle, KwakEngineCycle
from EngineArguments.default_arguments import get_default_kwargs
import pandas as pd
import numpy as np
import os


def get_kwak_data(EngineClass: KwakEngineCycle, design_args: dict, replication_mode: bool):
    default_kwargs = get_default_kwargs(EngineClass=EngineClass)
    total_kwargs = default_kwargs | {'_fuel_pump_pressure_factor_first_guess': 1.55,
                                     '_oxidizer_pump_pressure_factor_first_guess': 1.15,
                                     'is_frozen': False,
                                     'pressurant_heat_capacity_ratio': 1.667,
                                     'pressurant_molar_mass': 0.00399733779,
                                     'manual_oxidizer_density': 1126.1,
                                     'manual_fuel_density': 804.2,
                                     'replication_mode': replication_mode,
                                     }
    engine = EngineClass(**design_args, **total_kwargs)
    if issubclass(EngineClass, KwakFixGasGeneratorCycle):
        energy_source_list = [engine.turbine_propellant_mass, None]
    elif issubclass(EngineClass, KwakFixElectricPumpCycle):
        energy_source_list = [None, engine.battery.mass]
    return [
        engine.chamber_propellant_mass,
        *energy_source_list,
        engine.feed_system_mass,
        engine.tanks_mass,
        engine.pressurant.mass,
        engine.mass_kwak,
        engine.inverse_mass_ratio_kwak,
        engine.overall_specific_impulse,
        engine.ideal_delta_v_kwak,
    ]


design_args = {
    'thrust': 100e3,
    'combustion_chamber_pressure': 10e6,
    'burn_time': 300,
    'verbose': False
}

data = []
for EngineClass in [KwakFixElectricPumpCycle, KwakFixGasGeneratorCycle]:
    for replication_mode in [False, True]:
        engine_data = get_kwak_data(EngineClass=EngineClass, design_args=design_args, replication_mode=replication_mode)
        data.append(engine_data)
# for i in [0, 3]:
#     diff_list = [(1 - repl / kwak) * 100 if (kwak is not None and repl is not None) else None for kwak, repl in zip(data[i], data[i + 1]) ]
#     data.insert(i + 2, diff_list)
index = ['Kwak EP',
         'Repl. EP',
         'Kwak GG',
         'Repl. GG', ]
columns = ['CC Propellants [$kg$]',
           'GG Propellants [$kg$]',
           'Battery Pack [$kg$]',
           'Feed System [$kg$]',
           'Tanks [$kg$]',
           'Helium [$kg$]',
           'Total [$kg$]',
           'MR [-]',
           'Specific Impulse [$s$]',
           'Velocity Change [$m/s$]']

def round_to_significant_figures(x: float, n_figures: int):
    return round(x, n_figures - int(floor(log10(abs(x)))) - 1)

df = pd.DataFrame(data, index=index, columns=columns)
df = df.transpose()




def percentage_change(new_col, old_col):
    return (new_col / old_col - 1) * 100


for name in ['EP', 'GG']:
    df[f'Diff. {name}'] = percentage_change(df[f'Repl. {name}'], df[f'Kwak {name}'])

path = r'C:\Users\rvand\PycharmProjects\ThesisEP2\KwakFix\replication_table.xlsx'
with pd.ExcelWriter(path) as writer:
    df.to_excel(writer)
os.system(f'start EXCEL.EXE {path}')