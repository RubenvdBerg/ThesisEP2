from BaseEngineCycle.Cooling import CoolingChannels
from BaseEngineCycle.RP1Cooling import RP1CoolingChannels
from BaseEngineCycle.EngineCycle import EngineCycle
from BaseEngineCycle.SimpleCooling import Simple_CoolingChannels
import arguments as args

# for propellant_name in ('LH2_NASA', 'LO2_NASA', 'RP1_NASA'):
#     coolant = Coolant(propellant_name=propellant_name)
#     channels = CoolingChannels(coolant=coolant, total_heat_transfer=5E9, outlet_pressure=1E6, mass_flow=1000)
# channels = Simple_CoolingChannels(propellant_name='RP1_NASA', total_heat_transfer=11.68E6, outlet_pressure=10E6,
#                                   mass_flow=3.748969074623172 * 5.6)
# channels_rp1 = RP1CoolingChannels(total_heat_transfer=11.68E6, outlet_pressure=10E6, mass_flow=3.748969074623172 * 5.6)
# print(channels.outlet_temperature)
# print(channels_rp1.outlet_temperature)
# arguments = args.desgin_arguments.copy() | args.base_arguments.copy() | {'fuel_pump_specific_power': 15E3, 'oxidizer_pump_specific_power': 20E3}
# arguments['fuel_name'] = 'RP1_NASA'
# arguments['thrust'] = 10E3
# arguments['combustion_chamber_pressure'] = 10E6
# arguments['mass_mixture_ratio'] = 2.45
# arguments['is_frozen'] = False
# ec = EngineCycle(**arguments)
# print(ec)

# TDC1 kwargs

engine = EngineCycle(**args.tdc1_kwargs, **args.duel_pump_kwargs)
print(f'{engine.heat_exchanger.total_convective_heat_transfer*1e-6:.2f} MW', f'{engine.heat_exchanger.total_radiative_heat_transfer*1e-6:.2f} MW')
print(f'{engine.heat_exchanger.total_heat_transfer*1e-6:.2f} MW')
print(f'{engine.cooling_channels.increase_mass_specific_enthalpy*1e-6:.2f} MW/kg')
print(f'{engine.thrust_chamber.surface:.2f} m2')
engine.heat_exchanger.show_heat_flux()
engine.thrust_chamber.show_contour()

