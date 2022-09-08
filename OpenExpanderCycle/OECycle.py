import warnings
from dataclasses import dataclass

from BaseOpenCycle.OpenCycle import OpenEngineCycle


@dataclass
class OpenExpanderCycle(OpenEngineCycle):

    @property
    def turbine_mass_flow_initial_guess(self):
        return self.base_mass_flow * .03

    @property  # Override EngineCycle fuelflow -> increase pump requirements -> increase turbine mass flow -> iterate
    def main_fuel_flow(self):
        return 1 / (self.mass_mixture_ratio + 1) * self.chamber_mass_flow + self.turbine_mass_flow

    @property
    def turbine_inlet_temperature(self):
        if self.cooling_channel_section.outlet_temperature > self.maximum_wall_temperature:
            warnings.warn('Cooling outlet temperature cannot be higher than maximum wall temperature, cooling flow must be increased manually')
        return self.cooling_channel_section.outlet_temperature


