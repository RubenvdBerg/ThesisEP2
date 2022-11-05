from Verficiation.VV_test_heat_transfer import test_heat_transfer
from EngineCycles.BaseEngineCycle.HeatExchangerOMECA import RectangularOMECAHeatExchanger
import arguments as args
from math import pi
from dataclasses import dataclass, field
from functools import cached_property


# TODO: Review why it doesn't work

@dataclass
class LeonardiHeatExchanger(RectangularOMECAHeatExchanger):
    coolant_channel_throat_aspect_ratio: float = 3
    coolant_channel_height_input: float = field(init=False)

    @cached_property
    def section_coolant_channel_height(self):
        throat_circumference = self.thrust_chamber_section.get_radius(0) * 2 * pi
        throat_coolant_channel_width = throat_circumference / self.number_of_coolant_channels - self.coolant_channel_fin_thickness
        return throat_coolant_channel_width * self.coolant_channel_throat_aspect_ratio


# main_kwargs = {}
# for key, value in args.mira_kwargs.items():
#     if all(y not in key for y in ['turbine', 'exhaust', 'turbo']):
#         main_kwargs[key] = value

main_kwargs = {key: value for (key, value) in args.mira_kwargs.items() if
               all(y not in key for y in ['turbine', 'exhaust', 'turbo'])}
ht, plots = test_heat_transfer(engine_kwargs=main_kwargs,
                               number_of_coolant_channels=184,
                               chamber_wall_conductivity=365,
                               chamber_wall_thickness=1e-3,
                               heat_class=LeonardiHeatExchanger,
                               throat_area=0.116 ** 2 * pi / 4,
                               coolant_mass_flow=5.03,
                               coolant_inlet_temp=117,
                               coolant_inlet_pressure=76.8e5,
                               is_counter_flow=False,
                               verbose=True,
                               make_plots=True,
hot_gas_convective_heat_transfer_coefficient_mode='Bartz',
coolant_heat_transfer_coefficient_mode='SiederTate',
                               )

plots.plot_all()
