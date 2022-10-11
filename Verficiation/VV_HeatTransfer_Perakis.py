from Verficiation.VV_HeatTransfer_Denies import test_heat_transfer
import arguments as args
from math import radians, pi

main_kwargs = args.change_to_conical_nozzle(args.tcd1_kwargs, throat_half_angle=radians(25))
test_heat_transfer(engine_kwargs=main_kwargs,
                   number_of_coolant_channels=138,
                   chamber_wall_thickness=1e-3,
                   chamber_wall_conductivity=365,
                   coolant_mass_flow=2,
                   coolant_inlet_temp=30,
                   coolant_inlet_pressure=150e5,
                   _initial_flow_speed=10,
                   verbose=True,
                   counter_flow=True)


