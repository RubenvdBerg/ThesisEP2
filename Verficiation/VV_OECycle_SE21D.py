import arguments as args
from OpenExpanderCycle.OECycle import OpenExpanderCycle

engine = OpenExpanderCycle(**args.change_to_conical_nozzle(args.se_21d_kwargs), iterate=False)
data = [('Name', '[Unit]', 1, 1),
        ('Chamber Diameter', '[m]', engine.combustion_chamber.radius * 2, 0.985),
        ('Chamber Volume', '[m3]', engine.combustion_chamber.volume_incl_convergent, 1.029),
        ('Subs. Length', '[m]', engine.combustion_chamber.length + engine.nozzle.conv_length, 1.582),
        ('Throat Radius', '[m]', engine.nozzle.throat_radius, 0.286),
        ('Nozzle Length', '[m]', engine.nozzle.total_length, 2.548),
        # ('TC Length', '[m]', engine.thrust_chamber.length, 0),
        ('Turb. In. Temp.', '[K]', engine.cooling_channels.outlet_temperature, 506),
        ('Heat Transfer', '[MW]', engine.heat_transfer_section.total_heat_transfer*1e-6, 90),
        ('Turb. Mass Flow', '[kg/s]', engine.turbine_mass_flow, 10.2),
        ('Cool. Mass Flow', '[kg/s]', engine.cooling_flow, 16.394),
        # ('', '[]', 1, 1),
        # ('', '[]', 1, 1),
        # ('', '[]', 1, 1),
        ]

for name, unit, value, expected in data:
    print(f'{name:<18} \t {unit:<7} \t {value:<8.3f} \t {expected:.3f}')
