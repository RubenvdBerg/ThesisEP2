from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from EngineArguments.default_arguments import get_default_kwargs
from EngineFunctions.BaseFunctions import get_unit, get_symbol
from EngineFunctions.PlottingFunctions import adjust_values_to_prefix
from numpy import linspace
from typing import Iterable
import pandas as pd
import os

default_classes = (ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle)


def get_engines(engine_classes: Iterable[EngineCycle] = default_classes, **engine_kwargs):
    for EngineClass in engine_classes:
        default_kwargs = get_default_kwargs(EngineClass)
        total_kwargs = default_kwargs | engine_kwargs
        yield EngineClass(**total_kwargs)


def get_engines_mass_data(*engines: EngineCycle):
    for engine in engines:
        yield pd.DataFrame.from_dict(data=engine.combined_info, orient='index', columns=[type(engine).__name__])


def get_comparison_df(engine_classes: Iterable[EngineCycle] = default_classes, **engine_kwargs):
    return pd.concat(get_engines_mass_data(*get_engines(engine_classes, **engine_kwargs)), axis=1)


def get_comparison_dfs(input_variable: str, input_range: Iterable[float], input_nums: int, **kwargs):
    input_values = linspace(*input_range, input_nums)
    for input_value in input_values:
        kwargs[input_variable] = input_value
        yield get_comparison_df(**kwargs)


def get_comparison_excel(path: str, input_attribute: str, input_range: Iterable[float], input_nums: int,
                         input_prefix: str = '', **kwargs):
    input_values = linspace(*input_range, input_nums)
    unit = get_unit(input_attribute)
    symbol = get_symbol(input_attribute)
    with pd.ExcelWriter(path, engine='xlsxwriter') as writer:
        for input_value in input_values:
            kwargs[input_attribute] = input_value
            input_val = adjust_values_to_prefix([input_value], input_prefix)[0]
            sheet_name = f'{symbol}={input_val:.1f} {input_prefix}{unit}'
            get_comparison_df(**kwargs).to_excel(writer, sheet_name=sheet_name)
            workbook = writer.book
            worksheet = writer.sheets[sheet_name]
            format1 = workbook.add_format({'num_format': '0.0'})
            worksheet.set_column(1, 3, None, format1)
    os.system(f'start EXCEL.EXE {path}')



if __name__ == '__main__':
    kwargs1 = {
        'thrust': 100e3,
        'burn_time': 1200,
        # 'combustion_chamber_pressure': 10e6,
        'exit_pressure_forced': 0.002e6,
        # 'mass_mixture_ratio': 2.45,
        'expansion_ratio_end_cooling': 10,
        'maximum_wall_temperature': 900,
    }
    get_comparison_excel('test_excel.xlsx',input_attribute='combustion_chamber_pressure', input_range=(3e6, 10e6), input_nums=5, input_prefix='M',**kwargs1)
