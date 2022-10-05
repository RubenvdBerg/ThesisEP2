from Verficiation.VV_HeatTransfer_Denies import get_test_heat_transfer
import arguments as args
from math import radians, pi

get_test_heat_transfer(engine_kwargs=args.hyprob_kwargs,
                       throat_area=0.002879753,
                       number_of_coolant_channels=96,
                       chamber_wall_thickness=.9e-3,
                       chamber_wall_conductivity=365,
                       coolant_mass_flow=1.92,
                       coolant_inlet_temp=112.384,
                       coolant_inlet_pressure=15583600.,
                       verbose=True,
                       counter_flow=True)

channel_height = 2.45716897722*1e-3
channel_width = 4.53032220526*1e-3
channel_area = channel_width * channel_height
channel_equivalent_diameter = 2 * (channel_area / pi)**.5

Dh = 0.003186198562223674