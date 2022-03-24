import numpy as np
import pandas as pd
import pygmo as pg
from base_engine_cycle import EngineCycle
from base_gg_cycle import GasGeneratorCycle
from base_ep_cycle import ElectricPumpCycle
import matplotlib.pyplot as plt
import arguments as args


class EngineCycleOptimization:
    def __init__(self, thrust: float, burn_time: float, is_frozen: bool = False, verbose: bool = False,
                 cycle: type(EngineCycle) = ElectricPumpCycle):
        base_args = args.base_arguments
        if cycle == ElectricPumpCycle:
            extra_arguments = args.ep_arguments
        elif cycle == GasGeneratorCycle:
            extra_arguments = args.gg_arguments
        else:
            ValueError('cycle must be one of [ElectricPumpCycle, GasGeneratorCycle]')
        init_p_cc = 1E6
        self.ec = cycle(combustion_chamber_pressure=init_p_cc, thrust=thrust, burn_time=burn_time, is_frozen=is_frozen,
                        verbose=verbose, **base_args, **extra_arguments)

    def get_bounds(self):
        return [1.5, 1E5], [4, 3E7]


class InitialMassOpt(EngineCycleOptimization):
    def fitness(self, x):
        self.ec.mmr = x[0]
        self.ec.p_cc = x[1]
        self.ec.reiterate()
        return [self.ec.mass]


class DeltaVOpt(InitialMassOpt):
    def fitness(self, x):
        self.ec.mmr = x[0]
        self.ec.p_cc = x[1]
        self.ec.reiterate()
        return [-self.ec.ideal_delta_v]


def optimal_cycle_variables(thrust: float, burn_time: float, is_frozen: bool = False, verbose: bool = False,
                            cycle: type(EngineCycle) = ElectricPumpCycle, method: str = 'Nelder-Mead', n_pop: int = 1,
                            optimize_class: type(EngineCycleOptimization) = InitialMassOpt):
    prob = pg.problem(optimize_class(thrust, burn_time, is_frozen, verbose, cycle))
    algo = pg.algorithm(pg.scipy_optimize(method=method))
    pop = pg.population(prob, n_pop)
    pop_e = algo.evolve(pop)
    del prob
    return pop_e.champion_x, pop_e.champion_f


def plot_optimal_mmr_pcc(thrust: float, verbose: bool = False, method: str = 'Nelder-Mead', n_pop: int = 1,
                         burn_range: type(np.ndarray) = np.linspace(300, 1200, 15),
                         cycles: tuple = (GasGeneratorCycle, ElectricPumpCycle), mode: str = 'both'):
    data = []
    mode = mode.lower()
    assert mode in ['both', 'frozen', 'shifting']
    if mode == 'both':
        modes = [True, False]
    elif mode == 'frozen':
        modes = [True]
    elif mode == 'shifting':
        modes = [False]

    for is_frozen in modes:
        for cycle in cycles:
            mmr_data = []
            p_cc_data = []
            for burn_time in burn_range:
                prob = pg.problem(InitialMassOpt(thrust, burn_time, is_frozen, verbose, cycle))
                algo = pg.algorithm(pg.scipy_optimize(method=method))
                pop = pg.population(prob, n_pop)
                pop = algo.evolve(pop)
                mmr_opt, p_cc_opt = pop.champion_x
                mmr_data.append(mmr_opt)
                p_cc_data.append(p_cc_opt * 1E-6)
            mode_name = 'Frozen' if is_frozen else 'Shifting'
            cycle_name, color = ('EP', "blue") if cycle == ElectricPumpCycle else ('GG', 'red')
            linestyle = 'solid' if is_frozen else 'dashed'
            label = f'{cycle_name}_{mode_name.lower()}'
            data.append([p_cc_data, mmr_data, mode_name, cycle_name, color, linestyle, label])
            print(f'data made for {label}')
    fig, (ax1, ax2) = plt.subplots(2)
    fig.suptitle(f'Optimal MMR and Chamber Pressure \n{thrust * 1e-3} kN')
    ax1.set_ylabel('$p_{cc}$ [MPa]')
    ax2.set_ylabel('$MMR$ [-]')
    # ax1.set_ylim(top=9, bottom=6)
    plt.xlabel('$t_b$ [s]')
    for p_cc_data, mmr_data, mode, cycle_name, color, linestyle, label in data:
        ax1.plot(burn_range, p_cc_data, label=label, color=color, linestyle=linestyle)
        ax2.plot(burn_range, mmr_data, label=label, color=color, linestyle=linestyle)
    ax1.legend(ncol=2, loc=0)
    plt.show()
    plt.savefig('data/opt_vals_kwak1.png', dpi=600)




if __name__ == '__main__':
    # plot_optimal_mmr_pcc(thrust=10e3, mode='shifting', cycles=(ElectricPumpCycle,))
    # print(optimal_cycle_variables(thrust=10e3, burn_time=300, cycle=ElectricPumpCycle,
    #                               optimize_class=DeltaVOpt, verbose=True))
    print(optimal_cycle_variables(thrust=10e3, burn_time=300, cycle=GasGeneratorCycle,
                                  optimize_class=DeltaVOpt, verbose=True))
# m_f_cc: 0.8747795405015368
# m_fp: 0.8859693211184605, m_f_cool_act:0.0, m_f_cool_req: 0.011189780616923673
# m_fp: 0.886043371033555, m_f_cool_act:0.011189780616923661, m_f_cool_req: 0.011263830532018146
# m_fp: 0.8860438610690716, m_f_cool_act:0.011263830532018182, m_f_cool_req: 0.01126432056753477
# m_fp: 0.8860438643119491, m_f_cool_act:0.011264320567534791, m_f_cool_req: 0.011264323810412324
# (array([3.91515649e+00, 1.21205008e+06]), array([-10904.82098829]))
# GG:0.000000kg/s
# TU1:0.012816kg/s
# GG:0.012816kg/s
# TU1:0.012899kg/s
# Mass Flow Set
# (array([3.07357568e+00, 9.34322446e+06]), array([-12983.34290515]))
# (array([3.04705452e+00, 8.00392303e+06]), array([-11507.39899309]))
# (array([3.17706399e+00, 1.12213087e+07]), array([-12208.73827668]))