import matplotlib.pyplot as plt
import numpy as np
from rocketcea.cea_obj import CEA_Obj

# Pressures in MPa
cea_obj = CEA_Obj(fuelName='RP1_NASA', oxName='LO2_NASA')
pc_exit = .002
mrs = np.linspace(1.5, 4.0, 1000)

marker_switcher = {1: 'o', 5: 's', 10: '^', 20: 'x'}

def reduce_list(full_list, reduc: int = 50):
    return [x for i,x in enumerate(full_list) if i % reduc == 0]

for pc in list(marker_switcher)[::-1]:
    mpa_to_psia = 1.4503773800721814532099090803672E2
    pc_psia = round(pc * mpa_to_psia)
    eps = pc / pc_exit
    isps = [cea_obj.get_Isp(Pc=pc_psia, MR=mr, eps=eps, frozen=1, frozenAtThroat=1) for mr in mrs]
    max_isp = max(isps)
    max_mr = mrs[isps.index(max_isp)]
    reduced_mrs, reduced_isps = reduce_list(mrs), reduce_list(isps)
    data = [(mr, isp) for (mr, isp) in zip(reduced_mrs, reduced_isps) if isp > 50]
    data = list(zip(*data))

    plt.xlabel('Oxidizer-to-fuel ratio [-]')
    plt.ylabel('Chamber Specific Impulse [s]')
    plt.plot(*data, label=r'$p_{cc}$: ' + f'{pc} MPa', marker=marker_switcher[pc], color='black')
    plt.plot(max_mr, max_isp, linestyle='', color='r', marker='*')
    plt.text(max_mr-.15, max_isp + 2, r'O/F$_{opt}$:' + f'{max_mr:.2f}', color='r')
plt.ylim((310,385))
plt.legend()
plt.savefig('Isp_vs_MMR', dpi=900)
plt.show()