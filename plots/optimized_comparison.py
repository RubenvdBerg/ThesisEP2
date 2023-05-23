import matplotlib.pyplot as plt
import numpy as np

from new_new_optimization import optimize_engine
from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from EngineArguments.default_arguments import get_default_kwargs
from typing import Container
from EngineFunctions.PlottingFunctions import make_axis_string, get_class_color_marker, adjust_values_to_prefix
from EngineFunctions.BaseFunctions import get_symbol


def attribute_optimized_comparison(input_attribute: str,
                                   input_vals: Container[float],
                                   opt_attr: str = 'change_in_velocity',
                                   engine_classes: tuple = (ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle),
                                   input_prefix: str = '',
                                   opt_prefix: str = '',
                                   **engine_kwargs,
                                   ):
    acronym_switcher = {ElectricPumpCycle: 'EP', GasGeneratorCycle: 'GG', OpenExpanderCycle: 'OE'}

    def make_output(input_val: float):
        total_engine_kwargs = engine_kwargs | {input_attribute: input_val}
        output = {}
        for EngineClass in engine_classes:
            # Chamber pressure p_cc in MPa already
            p_cc, mmr, z = optimize_engine(
                EngineClass,
                attribute=opt_attr,
                total_kwargs=total_engine_kwargs | get_default_kwargs(EngineClass, False)
            )
            acronym = acronym_switcher[EngineClass]
            output[f'{acronym}_pcc'], output[f'{acronym}_mmr'], output[f'{acronym}_{opt_attr}'] = p_cc, mmr, z
            output[f'EP/{acronym}_{opt_attr}'] = output[f'EP_{opt_attr}'] / z
        return output

    outputs = [make_output(input_val) for input_val in input_vals]
    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    ax.set_xlabel(make_axis_string(input_attribute, input_prefix))
    output_symbol = get_symbol(opt_attr)
    ax.set_ylabel(make_axis_string(opt_attr, opt_prefix))
    ax2.set_ylabel('Optimal Chamber Pressure [MPa]')

    custom_legend = {'lines': [], 'labels':[]}
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
        adjust_values_to_prefix(output_vals, opt_prefix)
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
    plt.legend(handles=custom_legend['lines'], labels=custom_legend['labels'])
    plt.show()


if __name__ == '__main__':
    from plots.KwakPlots.Results_Comparison_RP1 import engine_kwargs

    attribute_optimized_comparison(input_attribute='thrust',
                                   input_vals=np.linspace(30e3, 100e3, 5),
                                   input_prefix='k',
                                   opt_prefix='k',
                                   engine_classes=[ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle],
                                   **engine_kwargs,
                                   )
