import json
from time import strftime
from typing import Container, Optional
from math import e
import numpy as np

from EngineArguments.default_arguments import get_default_kwargs
from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from new_new_optimization import optimize_engine
from optimized_comparison import plot_adj_attr_opt_compare_data


def get_sakama_electric_motor_specific_power(years: int) -> float:
    return 10 ** ((years - 2019) / 20) * 5.3E3


def get_tiede_battery_specific_energy(years: int, mode: int) -> float:
    mode_switcher = {
        1: lambda x: 1.922 * 10 ** (-28) * e ** (0.0344 * x),
        2: lambda x: 8.4692 * x - 16794,
        3: lambda x: 9.6056 * x - 19065,
    }
    func = mode_switcher[mode]
    return func(years)*3600


def make_adj_attr_opt_compare_data_years(engine_inputs: dict,
                                         input_attribute: str = 'years',
                                         input_vals: Container[float] = np.linspace(2000, 2500),
                                         opt_attr: str = 'change_in_velocity',
                                         engine_classes: tuple[EngineCycle] = (
                                                 ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle),
                                         input_filter: tuple[EngineCycle] = (),
                                         burn_times: tuple[float] = (300, 500, 1200),
                                         ignore_delta_p: bool = False,
                                         save_data: bool = True,
                                         log: bool = False,
                                         savepath: str = r"C:\\Users\\rvand\\PycharmProjects\\ThesisEP2\\data\\opt_compare\\"):
    if log:
        n = len(input_vals)
        start, stop = np.log10(input_vals[0]), np.log10(input_vals[-1])
        input_vals = np.logspace(start, stop, n)
    data_dict = {burn_time: {} for burn_time in burn_times}

    data_dict['info'] = {
        'input_attribute': input_attribute,
        'input_vals': list(input_vals),
        'opt_attr': opt_attr,
        'burn_times': burn_times,
        'engine_classes': [eng_class.__name__ for eng_class in engine_classes],
        'engine_inputs': engine_inputs,
        'log': log,
    }
    for burn_time in burn_times:
        engine_inputs['burn_time'] = burn_time

        def opt_return(EngineClass: EngineCycle, eng_kwargs: Optional[dict] = None):
            if eng_kwargs is None:
                eng_kwargs = engine_inputs
            if ignore_delta_p and EngineClass == ElectricPumpCycle:
                eng_kwargs['battery_specific_power'] = 1E20
            return optimize_engine(CycleClass=EngineClass,
                                   attribute=opt_attr,
                                   total_kwargs=get_default_kwargs(EngineClass, False) | eng_kwargs)

        if input_filter:
            filter_dict = {filter_engine_class: opt_return(filter_engine_class) for filter_engine_class in
                           input_filter}

        for EngineClass in engine_classes:
            outputs = {
                'chamber_pressure': [],
                'mixture_ratio': [],
                opt_attr: [],
            }
            data_dict[burn_time][EngineClass.__name__] = {}
            for input_val in input_vals:
                if EngineClass in input_filter:
                    output_vals = filter_dict[EngineClass]
                else:
                    input_dict = {
                        'battery_specific_energy': get_tiede_battery_specific_energy(input_val,1),
                        'electric_motor_specific_power': get_sakama_electric_motor_specific_power(input_val),
                    }
                    output_vals = opt_return(EngineClass, eng_kwargs=engine_inputs | input_dict)
                for val, key in zip(output_vals, outputs):
                    outputs[key].append(val)
            for key, values in outputs.items():
                data_dict[burn_time][EngineClass.__name__][key] = values

    if save_data:
        filename = strftime("%Y%m%d-%H%M%S") + '_DV_Opt_Comparison_YEARS'
        if ignore_delta_p:
            filename += '_Ignore_DP'
        with open(savepath + filename, 'w') as f:
            json.dump(data_dict, f, indent=6)
    return data_dict


if __name__ == '__main__':
    from plots.KwakPlots.Results_Comparison_RP1 import engine_kwargs
    # input_val = 2023 + 44.651680514950385
    # input_val = 2023 + 75.12506274555562
    # input_dict = {
    #     'battery_specific_energy': get_tiede_battery_specific_energy(input_val, 1),
    #     'electric_motor_specific_power': get_sakama_electric_motor_specific_power(input_val),
    # }

    data = make_adj_attr_opt_compare_data_years(
        input_attribute='years',
        input_vals=np.linspace(2067, 2098.1250627455556, 2),
        engine_classes=(ElectricPumpCycle,),
        # input_filter=(GasGeneratorCycle, OpenExpanderCycle),
        engine_inputs=engine_kwargs,
        burn_times=(300, 500, 1200),
        ignore_delta_p=True,
        log=False,
    )
    path = r'C:\Users\rvand\PycharmProjects\ThesisEP2\data\opt_compare\\'
    # filename = r'20230527-175835_DV_Opt_Comparison_YEARS_Ignore_DP'
    # with open(path + filename, 'r') as file:
    #     data = json.load(file)
    plot_adj_attr_opt_compare_data(
        data_dict=data,
        input_prefix='',
        output_prefix='k',
        # savefig='years',
    )
    # inputs = data['info']['input_vals']
    # dvs = {}
    # for burn_time in data:
    #     if burn_time != 'info':
    #         dvs[burn_time] = {}
    #         for engine_class, vals in data[burn_time].items():
    #             dvs[burn_time][engine_class] = vals['change_in_velocity']
    #
    # for burn_time in data:
    #     if burn_time != 'info':
    #         ep_dvs = dvs[burn_time]['ElectricPumpCycle']
    #         for conv_cycle in ['GasGeneratorCycle', 'OpenExpanderCycle']:
    #             lim = dvs[burn_time][conv_cycle][0]
    #             x = np.interp(lim, ep_dvs, inputs)
    #             # print(f'EP-cycle reaches the {conv_cycle} DV at {x} years for a burn time of {burn_time}')
    #             print(f'{conv_cycle}:{burn_time}:{x-2023}')