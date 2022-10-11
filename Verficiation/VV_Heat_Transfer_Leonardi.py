from Verficiation.VV_HeatTransfer_Denies import test_heat_transfer
from BaseEngineCycle.HeatTransferSection2 import DetailedHeatExchanger
import arguments as args
from math import pi
from dataclasses import dataclass, field
from functools import cached_property


@dataclass
class LeonardiHeatExchanger(DetailedHeatExchanger):
    coolant_channel_throat_aspect_ratio: float = 4
    coolant_channel_height_input: float = field(init=False)

    @cached_property
    def coolant_channel_height(self):
        throat_circumference = self.thrust_chamber.get_radius(0) * 2 * pi
        throat_coolant_channel_width = throat_circumference / self.number_of_coolant_channels - self.coolant_channel_fin_thickness
        return throat_coolant_channel_width * self.coolant_channel_throat_aspect_ratio


# main_kwargs = {}
# for key, value in args.mira_kwargs.items():
#     if all(y not in key for y in ['turbine', 'exhaust', 'turbo']):
#         main_kwargs[key] = value

main_kwargs = {key: value for (key, value) in args.mira_kwargs.items() if
               all(y not in key for y in ['turbine', 'exhaust', 'turbo'])}
test_heat_transfer(engine_kwargs=main_kwargs,
                   number_of_coolant_channels=160,
                   chamber_wall_conductivity=365,
                   chamber_wall_thickness=1e-3,
                   heat_class=LeonardiHeatExchanger,
                   throat_area=0.116 ** 2 * pi / 4,
                   coolant_mass_flow=5.03,
                   coolant_inlet_temp=117,
                   coolant_inlet_pressure=76.8e5,
                   counter_flow=False,
                   verbose=True,)
