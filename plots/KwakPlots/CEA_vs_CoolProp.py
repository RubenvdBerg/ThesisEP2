import pandas as pd
import matplotlib.pyplot as plt
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from EngineArguments.default_arguments import get_default_kwargs

from plots.KwakPlots.GG_vs_EP_Kwak import check_gg_vs_ep


def plot_cea_vs_coolprop_properties():
    base_path = r'C:\Users\rvand\PycharmProjects\ThesisEP2\data\Data\CEA_vs_CoolProp_'

    linestyle_switch = {
        650: '--',
        750: '-.',
        850: '-',
    }
    style_switch = {
        'CEA': ('o', 'b'),
        'CoolProp': ('s', 'r'),
    }

    title_switch = {
        'CP': ('Specific Heat Capacity', 'J/kg/K'),
        'MM': ('Molar Mass', 'g/mol'),
        'gamma': ('Heat Capacity Ratio', '-'),
    }

    for end_path in title_switch:
        path = base_path + end_path + '.csv'
        df = pd.read_csv(path, header=None, skiprows=(0, 1))
        fig, ax = plt.subplots()
        for column in df.iloc[:, :0:-1]:
            mode = 'CoolProp' if column < 4 else 'CEA'
            marker, color = style_switch[mode]
            temp = 650 + (column - 1) % 3 * 100
            linestyle = linestyle_switch[temp]
            label = f'{mode} - {temp} K'
            y = [x / 10 for x in df[0]]
            ax.plot(y, df[column],
                    marker=marker,
                    label=label,
                    color=color,
                    linestyle=linestyle)
        ax.set_xlabel(r'Pressure [MPa]')
        title, unit = title_switch[end_path]
        ax.set_ylabel(f'{title} [{unit}]')
        # Fix Custom Legend
        custom_lines = []
        custom_labels = []
        for temp in linestyle_switch:
            linestyle = linestyle_switch[temp]
            custom_lines.append(plt.Line2D([0], [0], linestyle=linestyle, color='black'))
            custom_labels.append(f'{temp} K')
        for mode in style_switch:
            marker, color = style_switch[mode]
            custom_lines.append(plt.Line2D([0], [0], color=color, marker=marker, linestyle="None"))
            custom_labels.append(mode)

        ncol = 2 if end_path == 'CP' else 1
        loc = (.1, .6) if end_path == 'CP' else 'best'
        ax.legend(custom_lines, custom_labels, ncol=ncol, loc=loc)
        plt.savefig(f'CEA_vs_CoolProp_{end_path}', dpi=1200)
        plt.show()


def check_performance_difference():
    from Results_Comparison_RP1 import engine_kwargs
    default_kwargs = get_default_kwargs(OpenExpanderCycle)
    engine_kwargs['turbine_maximum_temperature'] = 850
    engine_kwargs['combustion_chamber_pressure'] = 10e6
    total_kwargs = default_kwargs | engine_kwargs
    engine = OpenExpanderCycle(**total_kwargs)
    engine.print_masses()
    print(f'CEA {engine_kwargs["combustion_chamber_pressure"]/1e6:.0f} MPa Isp: {engine.overall_specific_impulse}')


if __name__ == '__main__':
    # plot_cea_vs_coolprop_properties()
    check_performance_difference()
    # check_gg_vs_ep()
