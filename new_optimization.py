import pygmo as pg

import scipy.optimize

from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle_DoublePump
from EngineArguments import arguments as args
from EngineArguments.default_arguments import get_default_kwargs


class EngineCycleOptimization:
    def __init__(self, cycle: type(EngineCycle), **kwargs):
        def_kwargs = get_default_kwargs(cycle, mass_mixture_ratio=False)
        self.kwargs = def_kwargs | kwargs
        self.cycle = cycle

    def get_bounds(self):
        if 'LH2' in self.kwargs['fuel_name']:
            return [4.5, 1E5], [8, 2E7]
        elif 'RP' in self.kwargs['fuel_name']:
            return [2, 1E5], [2.6, 1E7]


class DeltaVOpt(EngineCycleOptimization):
    def fitness(self, x):
        engine = self.cycle(mass_mixture_ratio=x[0], combustion_chamber_pressure=x[1], **self.kwargs)
        print(f'mmr: {x[0]:.2f}, pcc: {x[1]*1e-6:.1f}')
        return [-engine.ideal_delta_v]


class DeltaVOpt2(EngineCycleOptimization):
    def fitness(self, x):
        engine = self.cycle(mass_mixture_ratio=x[0], combustion_chamber_pressure=x[1], **self.kwargs)
        return [-engine.ideal_delta_v]

    def gradient(self, x):
        return pg.estimate_gradient_h(lambda y: self.fitness(y), x)


def optimal_cycle_variables(cycle: type(EngineCycle) = ElectricPumpCycle,
                            library: str = 'scipy',
                            n_pop: int = 1,
                            **kwargs):
    optimize_class = DeltaVOpt if library == 'scipy' else DeltaVOpt2
    prob = pg.problem(optimize_class(cycle, **kwargs))
    if library == 'nlopt':
        algo = pg.algorithm(pg.nlopt(solver='slsqp'))
    elif library == 'scipy':
        algo = pg.algorithm(pg.scipy_optimize(method='Nelder-Mead'))
    else:
        raise ValueError
    algo.set_verbosity(1)
    pop = pg.population(prob, n_pop)
    pop_e = algo.evolve(pop)
    del prob
    return pop_e.champion_x, pop_e.champion_f


if __name__ == '__main__':
    from plots.KwakPlots.Results_Comparison_RP1 import engine_kwargs
    design_variables = {
        'thrust': 100e3,
        'burn_time': 1200,
        'verbose': True,
    }
    engine_kwargs.pop('combustion_chamber_pressure')
    total_kwargs = engine_kwargs | design_variables
    print(optimal_cycle_variables(cycle=ElectricPumpCycle, **total_kwargs, library='scipy'))
