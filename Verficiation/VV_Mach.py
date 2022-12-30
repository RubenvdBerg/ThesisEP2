import matplotlib.pyplot as plt
from numpy import linspace
from EngineFunctions.IRTFunctions import get_local_mach, get_local_mach_nasa

num = 100
ar_range1 = list(linspace(1, 30, num))
ar_range2 = list(linspace(1, 2, num))
for is_subsonic, ar_range in zip([False, True], [ar_range1, ar_range2]):
    fig, ax = plt.subplots()
    for i in range(1, 5):
        y = 1 + i * .1
        if is_subsonic:
            print()
        for mach, style in zip([get_local_mach, get_local_mach_nasa], ['-', '-.']):
            m_range = [mach(local_area_ratio=ar,
                            is_subsonic=is_subsonic,
                            heat_capacity_ratio=y) for ar in ar_range]
            postfix = 'b4wind' if mach == get_local_mach else 'nasa'
            ax.plot(ar_range, m_range, label=f'y:{y} - {postfix}', linestyle=style)
    ax.set_xlabel('Area Ratio [-]')
    ax.set_ylabel('Mach [-]')
    plt.legend()
    plt.show()