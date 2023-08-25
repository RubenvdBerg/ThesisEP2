from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from plots.KwakPlots.isp_plot import get_data_dict


def graph_lambda_comparisons(attribute: str = 'mass_ratio', y_label: str = '$m_0/m_f$ [-]'):
    thrust = 100e3
    attribute_kwak = f'{attribute}_kwak'
    p_cc_range = tuple(range(6, 11))

    common_kwargs = {'thrusts': (thrust,),
                     'burn_times': (300, 390, 1200),
                     'names': ('EP',),
                     'p_cc_range': p_cc_range,
                     'exit_pressure_forced': 0.002e6}
    data_dict_kwak = get_data_dict(attributes=(attribute_kwak,),
                                   is_frozen=False,
                                   kwak=True,
                                   battery_fix_on=True,
                                   **common_kwargs)

    data_dict_repl = get_data_dict(attributes=(attribute_kwak,),
                                   is_frozen=False,
                                   kwak=True,
                                   battery_fix_on=False,
                                   **common_kwargs)

    data_dict_tool = get_data_dict(attributes=(attribute,),
                                   is_frozen=True,
                                   kwak=False,
                                   **common_kwargs)

    data_kwak = data_dict_kwak[attribute_kwak]['EP']
    data_repl = data_dict_repl[attribute_kwak]['EP']
    data_tool = data_dict_tool[attribute]['EP']

    fig, ax = plt.subplots()
    styles = ['-', '-.', '--']
    lines = []
    for color, data_dict in zip(['#ff9500', '#0C5DA5', '#07b547'], [data_tool, data_repl, data_kwak]):
        for style, (burn_time, data) in zip(styles, data_dict.items()):
            label = rf'{burn_time:.0f} s'
            line, = ax.plot(list(p_cc_range), data[thrust], linestyle=style, color=color,
                            label=label)
            lines.append(line)
    line1 = Line2D(p_cc_range, data[thrust])
    line1.set_linestyle('')
    ax.set_xlabel('$p_{cc}$ [MPa]')
    ax.set_ylabel(y_label)
    labels = [line.get_label() for line in lines]
    lines_ = [line1, *lines[:3], line1, *lines[3:6], line1, *lines[6:]]
    labels = [r"Final", *labels[:3], r"Init. Repl.", *labels[3:6], r'Adj. Repl.', *labels[6:]]
    ax.legend(lines_, labels, ncol=3, loc='lower right', bbox_to_anchor=(1.1, .8), shadow=True)
    path = rf'C:\Users\rvand\OneDrive\Documents\University of Technology Delft\Thesis\Figures\EP_Chapter_Graph_{attribute}.png'
    fig.savefig(path, dpi=1200)
    plt.show()


def make_comparison_tool_and_replication(data_dict_tool: dict, data_dict_repl: dict):
    p_cc_range = data_dict_repl['info']['p_cc_range']

    data_kwak = data_dict_repl[attribute_kwak]['EP']
    data_eta = data_dict_tool[attribute]['EP']

    fig, ax = plt.subplots()
    styles = ['-', '-.', '--']
    lines = []
    for color, data_dict in zip(['#ff9500', '#0C5DA5'], [data_eta, data_kwak]):
        for style, (burn_time, data) in zip(styles, data_dict.items()):
            label = rf'{burn_time:.0f} s'
            line, = ax.plot(list(p_cc_range), list(1 / x for x in data[thrust]), linestyle=style, color=color,
                            label=label)
            lines.append(line)
    line1 = Line2D(p_cc_range, data[thrust])
    line1.set_linestyle('')
    ax.set_xlabel('$p_{cc}$ [MPa]')
    ax.set_ylabel('$m_0/m_f$ [-]')
    labels = [line.get_label() for line in lines]
    lines_ = [line1, *lines[:3], line1, *lines[3:]]
    labels = [r"Final", *labels[:3], r"Init. Repl.", *labels[3:]]
    ax.legend(lines_, labels, ncol=2)
    fig.savefig(
        r'C:\Users\rvand\OneDrive\Documents\University of Technology Delft\Thesis\Figures\Diff_EP_Initial_Mass_Batt_Eta2.png')
    plt.show()


def lambda_comparison_kwak_vs_replication(data_dict_kwak, data_dict_repl):
    p_cc_range = data_dict_repl['info']['p_cc_range']
    data_kwak = data_dict_kwak[attribute]['EP']
    data_eta = data_dict_repl[attribute]['EP']

    fig, ax = plt.subplots()
    styles = ['-', '-.', '--']
    lines = []
    for color, data_dict in zip(['#ff9500', '#0C5DA5'], [data_eta, data_kwak]):
        for style, (burn_time, data) in zip(styles, data_dict.items()):
            label = rf'{burn_time:.0f} s'
            line, = ax.plot(list(p_cc_range), list(1 / x for x in data[thrust]), linestyle=style, color=color,
                            label=label)
            lines.append(line)
    line1 = Line2D(p_cc_range, data[thrust])
    line1.set_linestyle('')
    ax.set_xlabel('$p_{cc}$ [MPa]')
    ax.set_ylabel('$m_0/m_f$ [-]')
    labels = [line.get_label() for line in lines]
    lines_ = [line1, *lines[:3], line1, *lines[3:]]
    labels = [r"Init. Repl.", *labels[:3], r"Adj. Repl.", *labels[3:]]
    ax.legend(lines_, labels, ncol=2)
    fig.savefig(
        r'C:\Users\rvand\OneDrive\Documents\University of Technology Delft\Thesis\Figures\Diff_EP_Initial_Mass_Batt_Eta2.png')
    plt.show()


if __name__ == '__main__':
    graph_lambda_comparisons()
    graph_lambda_comparisons('ideal_delta_v', y_label='$\Delta V$ [m/s]')
