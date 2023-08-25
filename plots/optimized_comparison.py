import matplotlib.pyplot as plt
import numpy as np
from typing import Container, Optional
from time import strftime
import json

from optimization import optimize_engine
from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from EngineArguments.default_arguments import get_default_kwargs
from EngineFunctions.PlottingFunctions import make_axis_string, get_class_color_marker, adjust_values_to_prefix, \
    get_class_from_name
from EngineFunctions.BaseFunctions import get_symbol

acronym_switcher = {ElectricPumpCycle: 'EP', GasGeneratorCycle: 'GG', OpenExpanderCycle: 'OE'}


def attribute_optimized_comparison(input_attribute: str,
                                   input_vals: Container[float],
                                   opt_attr: str = 'change_in_velocity',
                                   engine_classes: tuple[EngineCycle] = (
                                           ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle),
                                   input_prefix: str = '',
                                   opt_prefix: str = '',
                                   savefig: Optional[str] = None,
                                   input_filter: tuple[EngineCycle] = (),
                                   **engine_kwargs,
                                   ):
    def opt_return(EngineClass: EngineCycle, eng_kwargs: Optional[dict] = None):
        if eng_kwargs is None:
            eng_kwargs = engine_kwargs
        return optimize_engine(CycleClass=EngineClass,
                               attribute=opt_attr,
                               total_kwargs=get_default_kwargs(EngineClass, False) | eng_kwargs)

    if input_filter:
        filter_dict = {filter_engine_class: opt_return(filter_engine_class) for filter_engine_class in input_filter}

    def make_output(input_val: float):
        output = {}
        for EngineClass in engine_classes:
            if EngineClass in input_filter:
                p_cc, mmr, z = filter_dict[EngineClass]
            else:
                # Chamber pressure p_cc in MPa already
                p_cc, mmr, z = opt_return(EngineClass, eng_kwargs=engine_kwargs | {input_attribute: input_val})
            acronym = acronym_switcher[EngineClass]
            output[f'{acronym}_pcc'], output[f'{acronym}_mmr'], output[f'{acronym}_{opt_attr}'] = p_cc, mmr, z
            output[f'EP/{acronym}_{opt_attr}'] = output[f'EP_{opt_attr}'] / z
        return output

    outputs = [make_output(input_val) for input_val in input_vals]
    input_vals = adjust_values_to_prefix(input_vals, si_prefix=input_prefix)
    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    ax.set_xlabel(make_axis_string(input_attribute, input_prefix))
    output_symbol = get_symbol(opt_attr)
    ax.set_ylabel(make_axis_string(opt_attr, opt_prefix))
    ax2.set_ylabel('Optimal Chamber Pressure [MPa]')

    custom_legend = {'lines': [], 'labels': []}
    for engine_class in engine_classes:
        acronym = acronym_switcher[engine_class]
        color, marker = get_class_color_marker(engine_class)
        pcc_vals = [output[f'{acronym}_pcc'] for output in outputs]
        ax2.plot(input_vals, pcc_vals,
                 label=r'$p_{cc}$' + f'{acronym}',
                 color=color,
                 marker=marker,
                 linestyle='--')
        output_vals = [output[f'{acronym}_{opt_attr}'] for output in outputs]
        output_vals = adjust_values_to_prefix(output_vals, opt_prefix)
        ax.plot(input_vals, output_vals,
                label=output_symbol + r'$_{' + f'{acronym}' + r'}$',
                color=color,
                marker=marker,
                linestyle='-')
        custom_legend['lines'].append(plt.Line2D([0], [0], color=color, marker=marker, linestyle="None"))
        custom_legend['labels'].append(f'{acronym}-cycle')

    custom_legend['lines'].append(plt.Line2D([0], [0], color='black', linestyle='--'))
    custom_legend['labels'].append(r'$p_{cc,opt}$')
    custom_legend['lines'].append(plt.Line2D([0], [0], color='black', linestyle='-'))
    custom_legend['labels'].append(output_symbol)
    ax.set_title(rf'$t_b$: {engine_kwargs["burn_time"]:.0f} s' + rf',$F_T$: {engine_kwargs["thrust"] * 1e-3:.0f} kN')
    plt.legend(handles=custom_legend['lines'], labels=custom_legend['labels'])
    if savefig:
        plt.savefig(f'{opt_attr}_opt_vs_{input_attribute}_{savefig}', dpi=1000)
    plt.show()


def make_adj_attr_opt_compare_data(input_attribute: str,
                                   input_vals: Container[float],
                                   engine_inputs: dict,
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
            filter_dict = {filter_engine_class: opt_return(filter_engine_class) for filter_engine_class in input_filter}

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
                    output_vals = opt_return(EngineClass, eng_kwargs=engine_inputs | {input_attribute: input_val})
                for val, key in zip(output_vals, outputs):
                    outputs[key].append(val)
            for key, values in outputs.items():
                data_dict[burn_time][EngineClass.__name__][key] = values

    if save_data:
        filename = strftime("%Y%m%d-%H%M%S") + '_DV_Opt_Comparison'
        if ignore_delta_p:
            filename += '_Ignore_DP'
        with open(savepath + filename, 'w') as f:
            json.dump(data_dict, f, indent=6)
    return data_dict


def plot_adj_attr_opt_compare_data(data_dict: dict,
                                   input_prefix: str = '',
                                   output_prefix: str = '',
                                   savefig: Optional[str] = None,
                                   log: Optional[bool] = None,
                                   output_attribute: Optional[str] = None,
                                   engine_classes: Optional[tuple[EngineCycle]] = None
                                   ):
    info = data_dict['info']
    input_attr = info['input_attribute']
    output_attr = info['opt_attr'] if output_attribute is None else output_attribute
    burn_times = [key for key in data_dict if key != 'info']
    engine_classes = info['engine_classes'] if engine_classes is None else [ec.__name__ for ec in engine_classes]
    input_vals_plot = adjust_values_to_prefix(info['input_vals'], si_prefix=input_prefix)
    if log is None:
        log = bool(info['log'])

    fig, ax = plt.subplots()
    if log:
        ax.set_xscale('log')
    xlabel = 'Year' if 'year' in input_attr else make_axis_string(input_attr, input_prefix)
    ylabel = make_axis_string(output_attr, output_prefix)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ax.set_ylabel('Optimal Chamber Pressure [MPa]')
    ax.set_title(rf'Optimized $\Delta V$ - $F_T$: {info["engine_inputs"]["thrust"] * 1e-3:.0f} kN')
    linestyles = ('dotted', 'dashed', 'solid', 'dashdot')

    # Plot lines
    for burn_time, linestyle in zip(burn_times, linestyles):
        for engine_class_name in engine_classes:
            EngineClass = get_class_from_name(engine_class_name)
            color, marker = get_class_color_marker(EngineClass)
            output_vals = data_dict[burn_time][engine_class_name][output_attr]
            output_vals_plot = adjust_values_to_prefix(output_vals, output_prefix)
            ax.plot(input_vals_plot, output_vals_plot,
                    color=color,
                    marker=marker,
                    linestyle=linestyle)

    # Make Custom Legend
    custom_legend = {'lines': [], 'labels': []}
    for engine_class_name in engine_classes:
        EngineClass = get_class_from_name(engine_class_name)
        color, marker = get_class_color_marker(EngineClass)
        acronym = acronym_switcher[EngineClass]
        custom_legend['lines'].append(plt.Line2D([0], [0], color=color, marker=marker, linestyle="None"))
        custom_legend['labels'].append(f'{acronym}-cycle')
    for burn_time, linestyle in zip(reversed(burn_times), reversed(linestyles[:-1])):
        custom_legend['lines'].append(plt.Line2D([0], [0], color='black', linestyle=linestyle))
        custom_legend['labels'].append(fr'$t_b$: {float(burn_time):.0f} s')
    plt.legend(handles=custom_legend['lines'], labels=custom_legend['labels'], ncol=2)

    if savefig:
        plt.savefig(f'{output_attr}_opt_vs_{input_attr}_{savefig}', dpi=1000)
    plt.show()


def adjusted_attribute_optimized_comparison(input_attr: str,
                                            input_vals: Container[float],
                                            engine_inputs: dict,
                                            output_attr: str = 'change_in_velocity',
                                            engine_classes: tuple[EngineCycle] = (
                                                    ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle),
                                            input_prefix: str = '',
                                            output_prefix: str = '',
                                            savefig: Optional[str] = None,
                                            input_filter: tuple[EngineCycle] = (),
                                            burn_times: tuple[float] = (300, 500, 1200),
                                            data_dict_path: Optional[str] = None,
                                            ignore_delta_p: bool = False,
                                            log: bool = False,
                                            output_attribute: Optional[str] = None,
                                            ):
    if data_dict_path is None:
        # Create Data
        data_dict = make_adj_attr_opt_compare_data(
            input_attribute=input_attr,
            input_vals=input_vals,
            opt_attr=output_attr,
            engine_classes=engine_classes,
            input_filter=input_filter,
            burn_times=burn_times,
            engine_inputs=engine_inputs,
            ignore_delta_p=ignore_delta_p,
            log=log,
        )
    else:
        with open(data_dict_path, 'r') as f:
            data_dict = json.load(f)

    plot_adj_attr_opt_compare_data(data_dict, input_prefix, output_prefix, savefig, log, output_attribute, engine_classes)


if __name__ == '__main__':
    from plots.KwakPlots.Results_Comparison_RP1 import engine_kwargs

    base_path = r'C:\Users\rvand\PycharmProjects\ThesisEP2\data\opt_compare\\'
    data_path = base_path + r'20230526-152429_DV_Opt_Comparison_Ignore_DP'
    # data_path = None
    adjusted_attribute_optimized_comparison(
        input_attr='battery_specific_energy',
        input_vals=np.linspace(150 * 3600, 100e6, 20),
        engine_classes=(ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle),
        input_filter=(GasGeneratorCycle, OpenExpanderCycle),
        engine_inputs=engine_kwargs,
        burn_times=(300, 500, 1200),
        input_prefix='',
        output_prefix='',
        savefig=None,
        data_dict_path=data_path,
        ignore_delta_p=True,
        log=True,
        output_attribute=None
    )
