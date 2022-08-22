import arguments as args
from OpenExpanderCycle.OECycle import OpenExpanderCycle

engine = OpenExpanderCycle(**args.le5a_kwargs_cnozzle, **args.open_arguments)
# engine.thrust_chamber.show_mach()
print(engine.cooling_channel_section.outlet_temperature)
# engine.thrust_chamber.show_contour()
# engine.heat_exchanger.show_heat_flux()
# print(engine.cooling_channels.outlet_temperature)
# print(engine.simple_specific_impulse)
# print(engine.thrust_chamber.length, engine.nozzle.exit_radius * 2)

