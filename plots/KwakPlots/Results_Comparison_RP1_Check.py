from isp_plot import get_engines
from EngineCycles.Abstract.EngineCycle import EngineCycle
from plots.Imaging.performance_image import make_performance_schematic


def get_engine(cycle: str, **kwargs):
    return list(get_engines(
        cycle,
        thrust=100e3,
        burn_time=1200,
        # p_cc_range=(10,),
        exit_pressure_forced=0.002e6,
        expansion_ratio_end_cooling=10,
        errors=True,
        verbose=True,
        is_frozen=True,
        **kwargs,
    ))[0]


def compare_engine_masses(*engines:EngineCycle):
    for attr in ['components_masses', 'aggregate_masses']:
        print('\n'+attr)
        componentss = [getattr(engine, attr) for engine in engines]
        combined_keys = list({key for components in componentss for key in components})
        combined_keys.sort()

        for key in combined_keys:
            masses = []
            for components in componentss:
                try:
                    masses.append(components[key])
                except KeyError:
                    masses.append(0)
            data = ",".join(f'{mass:>10.2f}' for mass in masses)
            print(f'{key:<25}: ' + data)

if __name__ == '__main__':
    oe = get_engine('ep', maximum_wall_temperature=900, p_cc_range=(3,),)
    oe2 = get_engine('ep', maximum_wall_temperature=900, p_cc_range=(4,),)
    oe3 = get_engine('ep', maximum_wall_temperature=900, p_cc_range=(5,),)
    oe4 = get_engine('ep', maximum_wall_temperature=900, p_cc_range=(6,),)
    oe5 = get_engine('ep', maximum_wall_temperature=900, p_cc_range=(7,),)
    oe6 = get_engine('ep', maximum_wall_temperature=900, p_cc_range=(8,),)
    oe7 = get_engine('ep', maximum_wall_temperature=900, p_cc_range=(9,),)
    oe8 = get_engine('ep', maximum_wall_temperature=900, p_cc_range=(10,),)
    # for engine in (oe, oe2,oe3,oe4,oe5,oe6,oe7,oe8):
    #     make_performance_schematic(engine)
    compare_engine_masses(oe, oe2,oe3,oe4,oe5,oe6,oe7,oe8)
    # gg = get_engine('gg')
    # ep = get_engine('ep')
    # compare_engine_masses(ep, oe , gg)
