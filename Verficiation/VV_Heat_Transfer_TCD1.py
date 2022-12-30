import numpy as np

from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineArguments import arguments as args
from numpy import linspace
import matplotlib.pyplot as plt
from math import radians

# file2 = pd.read_excel(filename +'.xlsx')


# TCD1 kwargs
main_kwargs = args.change_to_conical_nozzle(args.tcd1_kwargs, throat_half_angle=radians(25))
engine = EngineCycle(**main_kwargs, **args.duel_pump_kwargs)
main_kwargs['convective_coefficient_mode'] = 'Cornelisse'
engine_cornelis = EngineCycle(**main_kwargs, **args.duel_pump_kwargs)
print(engine.thrust_chamber.min_distance_from_throat)
# print(f'{engine.heat_exchanger.total_convective_heat_transfer*1e-6:.2f} MW', f'{engine.heat_exchanger.total_radiative_heat_transfer*1e-6:.2f} MW')
# print(f'{engine.heat_exchanger.total_heat_transfer*1e-6:.2f} MW')
# print(f'{engine.cooling_channels.increase_mass_specific_enthalpy*1e-6:.2f} MW/kg')
# print(f'{engine.thrust_chamber.surface:.2f} m2')
distance_tuple = engine.heat_transfer_section.distance_tuple
distances = list(linspace(*distance_tuple, 300))
heat_transfer_vals = [engine.heat_transfer_section.get_convective_heat_flux(distance) for distance in distances]
heat_transfer_vals_cornelis = [engine_cornelis.heat_transfer_section.get_convective_heat_flux(distance) for distance in
                               distances]
contour_vals = [engine.thrust_chamber.get_radius(distance) for distance in distances]

filename = r'C:\Users\rvand\PycharmProjects\ThesisEP2\BaseEngineCycle\Data\Perakis2021_HeatTransfer'
file1 = np.genfromtxt(fname=filename + '.txt', delimiter=',')
distances_exp = file1[1:, 0] * 1e-3 + distance_tuple[0]
values_exp = file1[1:, 1] * 1e6

fig, ax1 = plt.subplots()
ax1.set_xlabel('Distance from injection plate [mm]')
ticks3 = list(linspace(0 + distance_tuple[0], 1 + distance_tuple[0], 6))
ax1.set_xticks(ticks3)
ax1.set_xticklabels([f'{(x - distance_tuple[0]) * 1e3:.0f}' for x in ticks3])
ax1.set_xlim((0 + distance_tuple[0], 1 + distance_tuple[0]))
ax2 = ax1.twinx()

ax1.plot(distances, heat_transfer_vals, linestyle='-', label='Mod. Bartz', color='C0')
ax1.plot(distances, heat_transfer_vals_cornelis, linestyle=':', label='Cornelisse', color='C0')
ax1.plot(distances_exp, values_exp, linestyle='-', label='Perakis et al.', color='C1')
ax1.set_ylabel(r'Convective Heat Flux [$MW$/$m^2$]')
ticks1 = [x for x in range(0, 80000001, 10000000)]
ax1.set_yticks(ticks1)
ax1.set_yticklabels([f'{x * 1e-6:.1f}' for x in ticks1])
ax1.set_ylim((0., 80.e6))

ax2color = 'lightgrey'
ax2.plot(distances, contour_vals, linestyle='--', color=ax2color)
ax2.set_ylabel(r'Radius [mm]')
ticks2 = list(linspace(50e-3, 300e-3, 6))
ax2.set_yticks(ticks2)
ax2.set_yticklabels([f'{x * 1e3:.0f}' for x in ticks2])
ax2.set_ylim((50e-3, .3))
ax1.legend(loc='upper left')
plt.show()
engine.heat_transfer_section.show_heat_flux_coefficient()
engine.heat_transfer_section.show_heat_flux()
engine.heat_transfer_section.show_heat_transfer()
