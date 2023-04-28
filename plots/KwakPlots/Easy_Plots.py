from typing import Iterable, Optional, Literal
from operator import attrgetter
import matplotlib.pyplot as plt
from numpy import linspace

from EngineArguments.default_arguments import get_default_kwargs
from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineFunctions.BaseFunctions import get_unit, get_symbol, format_attr_name
from EngineFunctions.PlottingFunctions import make_axis_string, get_class_acronym, \
    get_class_color_marker, adjust_values_to_prefix, format_attr_name_for_legend, format_attr_name_for_axis_label
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle


def get_attribute_range(EngineClass: EngineCycle, input_attribute: str, input_range: tuple[float, float],
                        output_attributes: Iterable[str], num=5, **kwargs):
    default_kwargs = get_default_kwargs(EngineClass)
    total_kwargs = default_kwargs | kwargs
    input_values = linspace(*input_range, num)
    results = [input_values] + [tuple() for _ in range(len(output_attributes))]

    for input_value in input_values:
        total_kwargs = total_kwargs | {input_attribute: input_value}
        engine = EngineClass(**total_kwargs)
        for i, output_attribute in enumerate(output_attributes):
            results[i + 1] += attrgetter(output_attribute)(engine),
    return results


def make_engine_data_range(engine_classes: Iterable[EngineCycle], **kwargs):
    return {
        EngineClass: get_attribute_range(EngineClass, **kwargs) for EngineClass in engine_classes
    }


def double_input_data_range(input1_attribute: [str],
                            input2_attribute: [str],
                            input1_range: tuple[float, float],
                            input2_range: tuple[float, float],
                            num1: int = 5,
                            num2: int = 3, **kwargs):
    input2_values = linspace(*input2_range, num2)
    data_dict = {}
    for input2_val in input2_values:
        # label = f'{input2_attribute}:{input2_val}'
        kwargs[input2_attribute] = input2_val
        data = make_engine_data_range(input_attribute=input1_attribute,
                                      input_range=input1_range,
                                      num=num1,
                                      **kwargs, )
        data_dict[input2_val] = data
    return data_dict


def double_output_data_range(**kwargs):
    pass


def easy_plot(output_attribute: str,
              input_attribute: str = 'combustion_chamber_pressure',
              input_range: tuple[float] = (3e6, 10e6),
              engine_classes: Iterable[EngineCycle] = (ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle),
              num: int = 5, input_prefix: str = '', output_prefix: str = '',
              linestyle: str = '-',
              scale: bool = False,
              **engine_kwargs):
    fig, ax = plt.subplots()
    x_label = make_axis_string(input_attribute, input_prefix)
    y_label = make_axis_string(output_attribute, output_prefix)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    super_list = []

    for EngineClass in engine_classes:
        class_acronym = get_class_acronym(EngineClass)
        label = f'{class_acronym}-cycle'
        color, marker = get_class_color_marker(EngineClass)
        data = get_attribute_range(EngineClass, input_attribute, input_range, [output_attribute], num, **engine_kwargs)
        inputs, outputs = data
        input_vals = adjust_values_to_prefix(inputs, input_prefix)
        output_vals = adjust_values_to_prefix(outputs, output_prefix)
        if scale:
            super_list.append(output_vals)
            output_vals = [output_val / super_val for output_val, super_val in zip(output_vals, super_list[0])]
        ax.plot(input_vals, output_vals,
                color=color,
                linestyle=linestyle,
                marker=marker,
                label=label)
    ax.legend()
    plt.show()


def format_secondary_attribute_label_double_input(attribute: str, value: float, si_prefix: str):
    val = adjust_values_to_prefix([value], si_prefix)[0]
    unit = get_unit(attribute)
    try:
        symb = get_symbol(attribute)
    except KeyError:
        symb = attribute
    val_string = f'{val:.0f}{si_prefix}{unit}'
    return f'{symb}: {val_string}'


def double_input_plot(engine_classes: Iterable[EngineCycle],
                      input1_attribute: [str],
                      input2_attribute: [str],
                      input1_range: tuple[float, float],
                      input2_range: tuple[float, float],
                      output_attribute: str,
                      num1: int = 5,
                      num2: Literal[2, 3, 4] = 3,
                      input1_prefix: str = '',
                      input2_prefix: str = '',
                      output_prefix: str = '',
                      line_styles: Iterable[str] = ('-', '--', '-.', ':'),
                      x_log: bool = False,
                      y_log: bool = False,
                      savefig: Optional[str] = None,
                      scale: bool = False,
                      **engine_kwargs):
    # Get Data
    data_dict = double_input_data_range(input1_attribute=input1_attribute,
                                        input2_attribute=input2_attribute,
                                        input1_range=input1_range,
                                        input2_range=input2_range,
                                        num1=num1,
                                        num2=num2,
                                        engine_classes=engine_classes,
                                        output_attributes=[output_attribute],
                                        **engine_kwargs)

    # Set-up plot and plot lines
    fig, ax = plt.subplots()
    x_label = make_axis_string(input1_attribute, input1_prefix)
    y_label = make_axis_string(output_attribute, output_prefix)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if x_log:
        ax.set_xscale('log')
    if y_log:
        ax.set_yscale('log')

    for (input2_val, subdata_dict), linestyle in zip(data_dict.items(), line_styles):
        scale_list = []
        for EngineClass, data in subdata_dict.items():
            color, marker = get_class_color_marker(EngineClass)
            inputs, outputs = data

            inputs = adjust_values_to_prefix(inputs, input1_prefix)
            outputs = adjust_values_to_prefix(outputs, output_prefix)
            if scale:
                scale_list.append(outputs)
                outputs = [output / scale_list[0][0] for output in outputs]
            ax.plot(inputs, outputs, linestyle=linestyle, color=color, marker=marker)

    # Fix Custom Legend
    custom_lines = []
    custom_labels = []
    for EngineClass in engine_classes:
        color, marker = get_class_color_marker(EngineClass)
        acronym = get_class_acronym(EngineClass)
        custom_lines.append(plt.Line2D([0], [0], color=color, marker=marker, linestyle="None"))
        custom_labels.append(f'{acronym}-cycle')
    if num2 > 1:
        for input2_val, linestyle in zip(data_dict, line_styles):
            label = format_secondary_attribute_label_double_input(input2_attribute, input2_val, input2_prefix)
            custom_lines.append(plt.Line2D([0], [0], linestyle=linestyle, color='black'))
            custom_labels.append(label)
    ax.legend(custom_lines, custom_labels)

    if savefig:
        fig.savefig(savefig, dpi=1200)
    plt.show()
    return data_dict


def double_output_plot(engine_classes: Iterable[EngineCycle],
                       input_attribute: [str],
                       input_range: tuple[float, float],
                       output1_attribute: str,
                       output2_attribute: str,
                       num: int = 8,
                       input_prefix: str = '',
                       output1_prefix: str = '',
                       output2_prefix: str = '',
                       linestyles: tuple[str, str] = ('--', '-'),
                       twinx: bool = True,
                       legend_kwargs: Optional[dict] = None,
                       ylims: tuple = (None, None),
                       savefig: Optional[str] = None,
                       **engine_kwargs):
    x_label = make_axis_string(input_attribute, input_prefix)
    y1_label = make_axis_string(output1_attribute, output1_prefix)
    y2_label = make_axis_string(output2_attribute, output2_prefix)
    fig, ax = plt.subplots()
    ax2 = ax.twinx() if twinx else ax
    ax.set_xlabel(x_label)
    ax.set_ylabel(y1_label)
    ax2.set_ylabel(y2_label)
    ax.ticklabel_format(axis='y', style="sci", scilimits=(-2, 3))
    ax2.ticklabel_format(axis='y', style="sci", scilimits=(-2, 3))
    axes = [ax, ax2]

    for ylim, axis in zip(ylims, axes):
        if ylim is not None:
            axis.set_ylim(*ylim)

    output_attributes = [output1_attribute, output2_attribute]
    output_prefixes = [output1_prefix, output2_prefix]

    for EngineClass in engine_classes:
        inputs, *outputs_list = get_attribute_range(
            EngineClass, input_attribute, input_range, output_attributes, num, **engine_kwargs
        )
        input_vals = adjust_values_to_prefix(inputs, input_prefix)
        for axis, outputs, output_prefix, linestyle in zip(axes, outputs_list, output_prefixes, linestyles):
            class_acronym = get_class_acronym(EngineClass)
            label = f'{class_acronym}-cycle'
            color, marker = get_class_color_marker(EngineClass)
            output_vals = adjust_values_to_prefix(outputs, output_prefix)
            axis.plot(input_vals, output_vals,
                      color=color,
                      linestyle=linestyle,
                      marker=marker,
                      label=label, )

    # Fix Custom Legend
    custom_lines = []
    custom_labels = []
    for EngineClass in engine_classes:
        color, marker = get_class_color_marker(EngineClass)
        acronym = get_class_acronym(EngineClass)
        custom_lines.append(plt.Line2D([0], [0], color=color, marker=marker, linestyle="None"))
        custom_labels.append(f'{acronym}-cycle')
    for output_attribute, linestyle in zip([output1_attribute, output2_attribute], linestyles):
        label = format_attr_name_for_legend(output_attribute)
        custom_lines.append(plt.Line2D([0], [0], linestyle=linestyle, color='black'))
        custom_labels.append(label)
    if legend_kwargs is None:
        legend_kwargs = {}
    ax.legend(custom_lines, custom_labels, **legend_kwargs)

    if savefig:
        fig.savefig(savefig, dpi=1200)
    plt.show()


def compare_pressure(output_attribute: str,
                     engine_classes=[ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle],
                     input_attribute='combustion_chamber_pressure',
                     input_range=(3e6, 10e6),
                     input_prefix='M',
                     num=20,
                     output_prefix='',
                     **engine_kwargs):
    easy_plot(engine_classes=engine_classes,
              input_attribute=input_attribute,
              input_range=input_range,
              output_attribute=output_attribute,
              num=num,
              input_prefix=input_prefix,
              output_prefix=output_prefix,
              **engine_kwargs)


if __name__ == '__main__':
    from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
    from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
    from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
    from Results_Comparison_RP1 import engine_kwargs


    def turbine_explain(savefig: bool = False):
        double_output_plot(
            [GasGeneratorCycle, OpenExpanderCycle],
            'combustion_chamber_pressure',
            (3e6, 10e6),
            'turbine.inlet_flow_state.specific_heat_capacity',
            'turbine.inlet_flow_state.heat_capacity_ratio',
            input_prefix='M',
            legend_kwargs={'loc': 4, 'ncol': 2},
            num=8,
            savefig=r'Final_Plots\Turbine_Explain' if savefig else None,
            ylims=((3.6e3, 4.6e3), None),
            **engine_kwargs
        )
    def turbine_explain2(savefig: bool = False):
        double_output_plot(
            [GasGeneratorCycle, OpenExpanderCycle],
            'combustion_chamber_pressure',
            (3e6, 10e6),
            'turbine_effectivity',
            'turbine.mass_flow',
            input_prefix='M',
            legend_kwargs={'loc': 'best', 'ncol': 2},
            num=8,
            savefig=r'Final_Plots\Turbine_Explain2' if savefig else None,
            # ylims=((3.6e3, 4.6e3), None),
            **engine_kwargs
        )
    turbine_explain()
    # turbine_explain2()


    def feed_system_pressure(savefig: bool = False):
        double_output_plot(
            [ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle],
            'combustion_chamber_pressure',
            (3e6, 10e6),
            'total_thrust_chamber_mass',
            'feed_system_mass',
            input_prefix='M',
            legend_kwargs={'loc': 2, 'ncol': 2},
            savefig=r'Final_Plots\FS_TC_vs_Pressure' if savefig else None,
            num=8,
            ylims=(None, (0, 220)),
            # twinx=False,
            **engine_kwargs
        )


    def energy_source_pressure(savefig: bool = False):
        double_output_plot(
            [ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle],
            'combustion_chamber_pressure',
            (3e6, 10e6),
            'cc_prop_group_mass',
            'energy_source_mass',
            input_prefix='M',
            legend_kwargs={'loc': 'upper center', 'ncol': 2},
            savefig=r'Final_Plots\ES_vs_CCProp_Pressure' if savefig else None,
            num=5,
            ylims=((9600, 10300), (0, 700)),
            **engine_kwargs
        )


    def power_vs_anti_power_pressure():
        double_output_plot(
            [ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle],
            'combustion_chamber_pressure',
            (3e6, 10e6),
            'anti_power_mass',
            'power_mass',
            input_prefix='M',
            legend_kwargs={'loc': 9},
            # savefig=r'Final_Plots\Power_vs_Anti_Power_pressure',
            ylims=((9700, 10400), (100, 800)),
            **engine_kwargs
        )


    def power_vs_anti_power_burn_time():
        double_output_plot(
            [ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle],
            'burn_time',
            (300, 1200),
            'anti_power_ratio',
            'power_ratio',
            savefig=r'Final_Plots\Power_vs_Anti_Power_burn_time',
            legend_kwargs={'ncols': 2},
            ylims=((32.2, 33), None),
            **engine_kwargs
        )


    # turbine_explain(savefig=False)
    # energy_source_pressure(savefig=True)
    # feed_system_pressure(savefig=True)
    # power_vs_anti_power_burn_time()
    # power_vs_anti_power_pressure()

    def plot_pcc_exploration():
        easy_plot(
            output_attribute='anti_power_mass',
            input_prefix='M',
            **engine_kwargs
        )
        easy_plot(
            output_attribute='initial_mass',
            input_prefix='M',
            **engine_kwargs
        )
        easy_plot(
            output_attribute='power_mass',
            input_prefix='M',
            **engine_kwargs
        )


    def plot_tb_exploration():
        easy_plot(
            output_attribute='anti_power_ratio',
            input_attribute='burn_time',
            input_range=(300, 1200),
            input_prefix='',
            # scale=True,
            **engine_kwargs
        )
        easy_plot(
            output_attribute='initial_mass_ratio',
            input_attribute='burn_time',
            input_range=(300, 1200),
            input_prefix='',
            # scale=True,
            **engine_kwargs
        )
        easy_plot(
            output_attribute='power_ratio',
            input_attribute='burn_time',
            input_range=(300, 1200),
            input_prefix='',
            # scale=True,
            **engine_kwargs
        )


    # compare_pressure('initial_mass', **engine_kwargs)
    # compare_pressure('overall_specific_impulse', **engine_kwargs)
    # compare_pressure('chamber_propellant_mass', **engine_kwargs)
    # compare_pressure('tanks_mass', **engine_kwargs)
    # compare_pressure('feed_system_mass', **engine_kwargs)
    # compare_pressure('power_mass', **engine_kwargs)
    # compare_pressure('total_thrust_chamber_mass', **engine_kwargs)
    # compare_pressure('turbine.inlet_flow_state.specific_heat_capacity',
    #                  engine_classes=[GasGeneratorCycle, OpenExpanderCycle], **engine_kwargs, turbine_maximum_temperature=600)
    # compare_pressure('turbine.inlet_flow_state.heat_capacity_ratio',
    #                  engine_classes=[GasGeneratorCycle, OpenExpanderCycle], **engine_kwargs, turbine_maximum_temperature=600)
    # compare_pressure('turbine.inlet_pressure',engine_classes=[GasGeneratorCycle, OpenExpanderCycle], output_prefix='M', **engine_kwargs)
    # compare_pressure('turbine.inlet_flow_state.specific_heat_capacity',
    #                  engine_classes=[GasGeneratorCycle, OpenExpanderCycle], **engine_kwargs)
    # compare_pressure('turbine.inlet_flow_state.heat_capacity_ratio', engine_classes=[GasGeneratorCycle, OpenExpanderCycle],
    #                  **engine_kwargs)

    # compare_pressure('turbine.inlet_flow_state.molar_mass',engine_classes=[GasGeneratorCycle, OpenExpanderCycle], **engine_kwargs)
    # compare_pressure('fuel.mass', **engine_kwargs)
    # compare_pressure('oxidizer.mass', **engine_kwargs)
    # compare_pressure('propellant_related_mass', **engine_kwargs)
    #
    def make_turbine_inlet_comparison():
        double_input_plot([GasGeneratorCycle, OpenExpanderCycle],
                          input2_attribute='turbine_maximum_temperature',
                          input2_range=(850, 650),
                          num2=3,
                          input2_prefix='',
                          input1_attribute='combustion_chamber_pressure',
                          input1_range=(3e6, 10e6),
                          num1=5,
                          input1_prefix='M',
                          output_attribute='turbine.inlet_flow_state.molar_mass',
                          output_prefix='',
                          **engine_kwargs
                          )

    # double_input_plot([ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle],
    #                   input2_attribute='burn_time',
    #                   input2_range=(300, 1200),
    #                   num2=3,
    #                   input2_prefix='',
    #                   input1_attribute='combustion_chamber_pressure',
    #                   input1_range=(3e6, 10e6),
    #                   num1=5,
    #                   input1_prefix='M',
    #                   output_attribute='initial_mass',
    #                   output_prefix='',
    #                   **engine_kwargs
    #                   )
