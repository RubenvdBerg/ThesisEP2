import numpy as np
import pygmo as pg

import scipy.optimize

from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
import matplotlib.pyplot as plt
from EngineArguments import arguments as args


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
        self.ec.mass_mixture_ratio = x[0]
        self.ec.combustion_chamber_pressure = x[1]
        self.ec.update_cea()
        return [self.ec.mass]


class DeltaVOpt(EngineCycleOptimization):
    def fitness(self, x):
        self.ec.mass_mixture_ratio = x[0]
        self.ec.combustion_chamber_pressure = x[1]
        self.ec.update_cea()
        return [-self.ec.ideal_delta_v]


class DeltaVOpt2(EngineCycleOptimization):
    def fitness(self, x):
        self.ec.mass_mixture_ratio = x[0]
        self.ec.combustion_chamber_pressure = x[1]
        self.ec.update_cea()
        return [-self.ec.ideal_delta_v]

    def gradient(self, x):
        return pg.estimate_gradient_h(lambda y: self.fitness(y), x)


def optimal_cycle_variables_scipy(thrust: float, burn_time: float, is_frozen: bool = False, verbose: bool = False,
                                  cycle: type(EngineCycle) = ElectricPumpCycle, method: str = 'Nelder-Mead',
                                  n_pop: int = 1,
                                  optimize_class: type(EngineCycleOptimization) = InitialMassOpt):
    prob = pg.problem(optimize_class(thrust, burn_time, is_frozen, verbose, cycle))
    algo = pg.algorithm(pg.scipy_optimize(method=method))
    pop = pg.population(prob, n_pop)
    pop_e = algo.evolve(pop)
    del prob
    return pop_e.champion_x, pop_e.champion_f


def optimal_cycle_variables_nlopt(thrust: float, burn_time: float, cycle: type(EngineCycle), is_frozen: bool = False,
                                  verbose: bool = False, solver: str = 'slsqp', n_pop: int = 1,
                                  optimize_class: type(EngineCycleOptimization) = DeltaVOpt2):
    prob = pg.problem(optimize_class(thrust, burn_time, is_frozen, verbose, cycle))
    algo = pg.algorithm(pg.nlopt(solver=solver))
    algo.set_verbosity(1)
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


def set_attribute_cycle(attribute: str, cycle_type: str):
    cycle_types = ['ep', 'gg', 'ex', 'sc']
    if cycle_type == 'ep':
        cycle = ElectricPumpCycle
        extra_arguments = args.ep_arguments
    elif cycle_type == 'ex':
        raise NotImplementedError
    elif cycle_type == 'gg':
        cycle = GasGeneratorCycle
        extra_arguments = args.gg_arguments
    elif cycle_type == 'sc':
        raise NotImplementedError
    else:
        raise ValueError(f'Invalid cycle_type. Must be one of {cycle_types}')

    attribute_options = ['mass', 'dv']
    if attribute == 'mass':
        def attribute_function(x):
            return x.mass
    elif attribute == 'dv':
        def attribute_function(x):
            return -x.ideal_delta_v
    else:
        raise ValueError(f'Invalid attribute. Must be one of {attribute_options}')

    return cycle, extra_arguments, attribute_function


def fast_optimize(cycle_type: str, thrust: float, burn_time: float, is_frozen: bool = False,
                  verbose: bool = False, attribute: str = 'dv', x0=np.array([10, 25]), bounds=((1, 300), (15, 40)),
                  eps: float = 1e-2, basin_hop: int = 0, full_output: bool = False):
    base_args = args.base_arguments_o
    cycle, extra_arguments, attribute_function = set_attribute_cycle(attribute, cycle_type)

    base_args.pop('is_frozen')
    base_args['expansion_ratio'] = 50
    base_args['expansion_ratio_end_cooling'] = 5
    def optimize_function(x):
        engine_cycle = cycle(
            combustion_chamber_pressure=x[0] * 1e5,
            mass_mixture_ratio=x[1] * 1e-1,
            thrust=thrust,
            burn_time=burn_time,
            is_frozen=is_frozen,
            verbose=verbose,
            **base_args,
            **extra_arguments
        )
        return attribute_function(engine_cycle)

    if basin_hop:
        solution = scipy.optimize.basinhopping(
            func=optimize_function,
            x0=x0,
            minimizer_kwargs={
                'method': 'SLSQP',
                'bounds': bounds,
                'options': {'eps': eps}
            },
            niter=0,
            stepsize=2,
            disp=verbose,
            T=10
        )
    else:
        solution = scipy.optimize.minimize(fun=optimize_function,
                                           x0=x0,
                                           method='SLSQP',
                                           bounds=bounds,
                                           options={'eps': eps})
    if full_output:
        return solution
    return abs(solution.fun)


if __name__ == '__main__':
    solution1 = fast_optimize(cycle_type='gg',
                              thrust=10E3,
                              burn_time=600,
                              is_frozen=True,
                              verbose=False,
                              attribute='dv',
                              full_output=True
                              )
    print(solution1)
# plot_optimal_mmr_pcc(thrust=10e3, mode='shifting', cycles=(ElectricPumpCycle,))
# print(optimal_cycle_variables(thrust=10e3, burn_time=300, cycle=ElectricPumpCycle,
#                               optimize_class=DeltaVOpt, verbose=True))
# print(optimal_cycle_variables_scipy(thrust=10e3, burn_time=300, cycle=GasGeneratorCycle,
#                                     optimize_class=DeltaVOpt, verbose=False, is_frozen=False))
# for eps in [8e-3, 9e-3, 1.1e-2, 1.2e-2]:
#     time1 = perf_counter()
#     solution = fast_optimize(cycle=GasGeneratorCycle,
#                              thrust=10E3,
#                              burn_time=300,
#                              is_frozen=False,
#                              verbose=False,
#                              attribute='dv',
#                              eps=eps
#                              )
#     print(solution.fun)
#     print(perf_counter() - time1)
# print(optimal_cycle_variables_nlopt(thrust=10e3, burn_time=300, cycle=GasGeneratorCycle, n_pop=5,
#                                     optimize_class=DeltaVOpt2, verbose=False, is_frozen=False))
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

# GG 10 kN, 300 s
# Frozen (array([2.74903818e+00, 1.67141456e+07]), array([-10676.60988527]))
# Equilibrium (array([3.17406705e+00, 1.17634030e+07]), array([-12119.17445254]))
# (array([3.04101662e+00, 9.26083565e+06]), array([-12369.45333834]))
# (array([3.06947409e+00, 6.76381057e+06]), array([-12534.81489996]))
# (array([3.0454796e+00, 7.9300087e+06]), array([-11421.38748787]))

# (array([3.72961197e+00, 2.07122371e+07]), array([-11035.28802239]))
# lowest_optimization_result:      fun: -12609.126954853742
#      jac: array([0.00091584, 0.08859247])
#  message: 'Optimization terminated successfully'
#     nfev: 27
#      nit: 5
#     njev: 5
#   status: 0
#  success: True
#        x: array([43.71614877, 29.48556702])

# iter 0, stepsize 1e-3
#  lowest_optimization_result:      fun: -12609.126956025118
#      jac: array([0.00404969, 0.00293872])
#  message: 'Optimization terminated successfully'
#     nfev: 52
#      nit: 13
#     njev: 13
#   status: 0
#  success: True
#        x: array([43.72118336, 29.48551076])

# iter 0 stepsize 1e-2
# lowest_optimization_result:      fun: -12609.126783630098
#     jac: array([-0.00388071,  0.03883523])
# message: 'Optimization terminated successfully'
#    nfev: 24
#     nit: 6
#    njev: 6
#  status: 0
# success: True
#       x: array([43.69687819, 29.48204892])
