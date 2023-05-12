import pandas as pd
import matplotlib.pyplot as plt

base_path = r'C:\Users\rvand\PycharmProjects\ThesisEP2\data\Data\CEA_vs_CoolProp_'

end_paths = ['CP', 'MM', 'gamma']

linestyle_switch = {
    650: '-',
    750: '-.',
    850: '--',
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

for end_path in end_paths:
    path = base_path + end_path + '.csv'
    df = pd.read_csv(path, header=None, skiprows=(0, 1))
    fig, ax = plt.subplots()
    for column in df.iloc[:, 1:]:
        mode = 'CoolProp' if column < 4 else 'CEA'
        marker, color = style_switch[mode]
        temp = 650 + (column - 1) % 3 * 100
        linestyle = linestyle_switch[temp]
        label = f'{mode} - {temp} K'
        y = [x/10 for x in df[0]]
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
    for mode in style_switch:
        marker, color = style_switch[mode]
        custom_lines.append(plt.Line2D([0], [0], color=color, marker=marker, linestyle="None"))
        custom_labels.append(mode)
    for temp in linestyle_switch:
        linestyle = linestyle_switch[temp]
        custom_lines.append(plt.Line2D([0], [0], linestyle=linestyle, color='black'))
        custom_labels.append(f'{temp} K')
    ax.legend(custom_lines, custom_labels)
    plt.show()
