
from KwakFix.KwakFixCycles import KwakFixElectricPumpCycle
from isp_plot import plot_mr_kwak_ratio, get_data_dict
# import scienceplots
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
# import matplotlib as mpl
#
#
# plt.style.use('science')
# mpl.use('Qt5Agg')
# mpl.rcParams['axes.unicode_minus'] = False
thrust = 100e3
attribute = 'mass_ratio_kwak'
burn_times = 300, 390, 1200
data = [
    get_data_dict(attributes=(attribute,),
                  thrusts=(thrust,),
                  burn_times=burn_times,
                  names=('EP',),
                  kwak=True,
                  battery_fix_on=x,
                  p_cc_range=tuple(range(6,11)))
    for x in [True, False]
]
data_dict_kwak, data_dict_eta = data
p_cc_range = data_dict_kwak['info']['p_cc_range']

data_kwak = data_dict_kwak[attribute]['EP']
data_eta = data_dict_eta[attribute]['EP']

fig, ax = plt.subplots()
styles = ['-','-.','--']
lines = []
for color, data_dict in zip(['#ff9500', '#0C5DA5'],[data_eta, data_kwak]):
    for style, (burn_time,data) in zip(styles, data_dict.items()):
        label = rf'{burn_time:.0f} s'
        line, = ax.plot(list(p_cc_range), list(data[thrust]), linestyle=style, color=color, label=label)
        lines.append(line)
line1 = Line2D(p_cc_range, data[thrust])
line1.set_linestyle('')
ax.set_xlabel('$p_{cc}$ [MPa]')
ax.set_ylabel('Mass Ratio [-]')
labels = [line.get_label() for line in lines]
lines_ = [line1, *lines[:3], line1, *lines[3:]]
labels = [r"Incl. $\eta$'s", *labels[:3],  r"Excl. $\eta$'s", *labels[3:]]
ax.legend(lines_, labels, ncol=2)
fig.savefig(r'C:\Users\rvand\OneDrive\Documents\University of Technology Delft\Thesis\Figures\Diff_EP_Initial_Mass_Batt_Eta.png')
plt.show()
