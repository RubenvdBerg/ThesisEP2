from scipy.optimize import minimize
from EngineCycles.Abstract.EngineCycle import EngineCycle


def engine_opt_func(params, CycleClass: EngineCycle, attribute: str, total_kwargs: dict, is_max:bool = True):
    print(f'{params[0]:.3f} MPa, {params[1]:.3f}')
    total_kwargs['combustion_chamber_pressure'] = params[0] * 1e6
    total_kwargs['mass_mixture_ratio'] = params[1]
    multiplier = -1 if is_max else 1
    return multiplier * getattr(CycleClass(**total_kwargs), attribute)


def optimize_engine(CycleClass: EngineCycle, attribute: str, total_kwargs: dict, x0: tuple = (3, 2.4),
                    bounds: tuple = ((1, 10), (1.5, 3.5)), tol: float = 1e-2, is_max: bool = True):
    res = minimize(
        fun=lambda x: engine_opt_func(x, CycleClass, attribute, total_kwargs),
        x0=x0,
        bounds=bounds,
        method='Nelder-Mead',
        tol=tol,
    )
    fun_val = -res.fun if is_max else res.fun
    return *res.x, fun_val


if __name__ == '__main__':
    from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
    from EngineArguments.default_arguments import get_default_kwargs
    from plots.KwakPlots.Results_Comparison_RP1 import engine_kwargs

    total_kwargs = get_default_kwargs(ElectricPumpCycle) | engine_kwargs
    total_kwargs['burn_time'] = 1200
    print(
        optimize_engine(ElectricPumpCycle, 'ideal_delta_v', total_kwargs)
    )
