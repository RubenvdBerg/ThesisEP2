
import scienceplots
from rocketcea.cea_obj import CEA_Obj
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# plt.style.use('science')
# mpl.use('Qt5Agg')

cea = CEA_Obj(fuelName='RP1_NASA', oxName='LO2_NASA')


def get_t_comb(pc, MR):
    Pc = 14.503773800722 * pc
    tc_rankine = cea.get_Tcomb(Pc=Pc, MR=MR)
    return tc_rankine / 1.8


mmrs = [x * 1e-3 for x in range(250, 451, 15)]
for pc, marker in zip((100, 10, 1), ('^', 'o', 's')):
    temps = [get_t_comb(pc, mmr) for mmr in mmrs]
    plt.plot(mmrs, temps, label=f'CEA {pc} bar', color='#0C5DA5', marker=marker,
             markevery=range(1, len(mmrs) - 1))
temps2 = [1550.3 * mmr + 409.3 for mmr in mmrs]
plt.plot(mmrs, temps2, label='Choi et al.', linestyle='--', color='r')

data = ((922, .342, 'H-1'), (916, .33, 'RS-27'), (1062, .416, 'F-1'), (843.8, .297, 'S-4'), (900, .32, 'Kwak'))
# data = ((922, .342, 'H-1 (41.3 bar)'), (916, .33, 'RS-27 (48.7 bar)'), (1062, .416, 'F-1 (77.6 bar)'), (843.8, .297, 'S-4 (46 bar)'), (900, .32, 'Kwak (30-100 bar)'))

for temp, mmr, label in data:
    plt.plot(mmr, temp, linestyle='', marker='o', color='black')
    if 'H-1' in label:
        plt.text(mmr + .005, temp, label)
    else:
        plt.text(mmr + .005, temp - 20, label)

gg_pressures = {'F-1 ': '77.6', 'H-1 ': '41.3', 'RS-27': '48.7', 'Kwak': '30-100', 'S-4 ': '4.60'}
text = 'Engine: $p_{gg}$[bar]\n'
spaces = 8
for key, val in gg_pressures.items():
    if 'Kwak' in key or 'RS-27' in key:
        text += f'{key:<{spaces-2}}: {val}\n'
    else:
        text += f'{key:<{spaces}}: {val}\n'



plt.subplots_adjust(right=0.73)
plt.text(.75, .5, text, fontsize=14,transform=plt.gcf().transFigure)
plt.xlabel('O/F-ratio [-]', fontsize=14)
plt.ylabel('Combustion Temperature [K]', fontsize=14)
plt.legend()
plt.savefig('Choi2.png', dpi=1200)
plt.show()
