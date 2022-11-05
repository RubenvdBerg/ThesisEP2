import pygmo as pg
import numpy as np
from EngineCycles.GasGeneratorCycle.GGCycle import GasGeneratorCycle
import matplotlib.pyplot as plt


base_arguments = {
    'thrust': 100E3, 'burn_time': 300, 'combustion_chamber_pressure': 7.0E6, 'oxidizer_name': 'LO2_NASA',
    'fuel_name': 'RP1_NASA', 'is_frozen': True, 'exit_pressure': .002E6, 'max_acceleration': 4.5,
    'heat_ratio_pressurant': 1.667, 'mass_mixture_ratio': 4., 'pressurant_initial_pressure': 27E6,
    'pressurant_final_pressure': 5E6,
    'oxidizer_initial_pressure': .4E6, 'fuel_initial_pressure': .25E6, 'fuel_pump_pressure_factor': 1.55,
    'oxidizer_pump_pressure_factor': 1.15, 'pressurant_gas_constant': 2080,
    'pressurant_initial_temperature': 100, 'oxidizer_pump_efficiency': .66, 'fuel_pump_efficiency': .61,
    'pressurant_margin_factor': 1.1, 'pressurant_tank_structural_factor': 1.2,
    'propellant_margin_factor': 1.01, 'tanks_structural_factor': 2.5, 'ullage_volume_factor': 1.08,
    'oxidizer_density': 1126.1, 'fuel_density': 804.2, 'tanks_material_density': 2850,
    'pressurant_tank_material_density': 4430, 'tanks_yield_strength': 250E6,
    'pressurant_tank_yield_strength': 1100E6, 'verbose': True, 'kwak_fix': True
}

gg_arguments = {
    'gg_gas_specific_heat': 2024.7, 'heat_ratio_gg_gas': 1.16, 'mass_mixture_ratio_gg': 0.320,
    'turbine_pressure_ratio': 27, 'gas_constant_gg_gas': 274.1, 'turbine_inlet_temperature': 900,
    'gas_generator_stay_time': 10E-3, 'turbopump_specific_power': 13.5E3, 'turbine_efficiency': .52,
    'gg_structural_factor': 2.5, 'gg_material_density': 8220, 'gg_yield_strength': 550E6, 'gg_thrust_contribution': 0.01

}

class GGCycleOptimization:
    def __init__(self, arguments=base_arguments):
        self.ec = GasGeneratorCycle(**arguments, **gg_arguments)

    def fitness(self, x):
        self.ec.combustion_chamber_pressure = x[0]
        return [self.ec.mass]

    def get_bounds(self):
        return [1E6], [1E8]


def optimal_p_cc_gg(arguments, method, n_pop):
    prob = pg.problem(GGCycleOptimization(arguments=arguments))
    algo = pg.algorithm(pg.scipy_optimize(method=method))
    pop = pg.population(prob, n_pop)
    pop = algo.evolve(pop)
    return pop.champion_x[0], pop.champion_f[0]


if __name__ == '__main__':
    pressure_range = np.logspace(6, 9, 20)
    data = []
    for pcc in pressure_range:
        base_arguments['combustion_chamber_pressure'] = pcc
        gg = GasGeneratorCycle(**base_arguments, **gg_arguments)
        data.append(gg.mass)
    ft, mmr, tb = base_arguments["thrust"], base_arguments["mass_mixture_ratio"], base_arguments["burn_time"]
    frozen = 'FrozenAtThroat' if base_arguments["is_frozen"] else "Shifting Equilibrium"
    plt.plot(pressure_range, data)
    plt.axvline(pressure_range[data.index(min(data))], color='red')
    plt.title(r'GG ' + frozen + '\n'+f'MMR={mmr:.2f},' + r' $t_b$='+f'{tb:.0f} s')
    plt.ylabel(r'Initial Mass [kg]')
    plt.xlabel(r'Chamber Pressure [Pa]')
    plt.xscale('log')
    plt.show()

