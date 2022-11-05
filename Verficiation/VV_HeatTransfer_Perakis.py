import matplotlib.pyplot as plt

from Verficiation.VV_test_heat_transfer import test_heat_transfer
from EngineCycles.BaseEngineCycle.HeatExchanger import HeatExchanger
import arguments as args
from math import radians
import numpy as np



main_kwargs = args.change_to_conical_nozzle(args.tcd1_kwargs, throat_half_angle=radians(25))
heattransfer, plots = test_heat_transfer(engine_kwargs=main_kwargs,
                                         make_plots=True,
                                         number_of_coolant_channels=138,
                                         heat_class=HeatExchanger,
                                         maximum_wall_temp=800,
                                         coolant_mass_flow=23.14758424727635,
                                         coolant_inlet_temp=30,
                                         coolant_inlet_pressure=150e5,
                                         _initial_flow_speed=10,
                                         verbose=True,
                                         is_counter_flow=False, )

def extra_func(ax: plt.Axes, data:dict):
    min_distance = heattransfer.thrust_chamber_section.min_distance_from_throat
    filename = r'../data/Data\Perakis2021_HeatTransfer'
    file1 = np.genfromtxt(fname=filename + '.txt', delimiter=',')
    distances_exp = file1[1:, 0] * 1e-3 + min_distance
    values_exp = file1[1:, 1] * 1e6
    ax.plot(distances_exp, values_exp, label='q_Perakis', color='pink', linestyle='--')

plots.plot_flux(extra_func=extra_func)
plots.plot_coolant()