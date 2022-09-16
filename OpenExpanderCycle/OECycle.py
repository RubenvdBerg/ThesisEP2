import warnings
from dataclasses import dataclass

from BaseOpenCycle.OpenCycle import OpenEngineCycle


@dataclass
class OpenExpanderCycle(OpenEngineCycle):

    @property
    def default_turbine_flow_check_state(self):
        return self.cooling_channel_section.outlet_flow_state

    @property
    def turbine_mass_flow_initial_guess(self):
        return self.base_mass_flow * .03

    @property  # Overrides flow exiting fuel tank -> increases pump requirements -> increases turbine mass flow -> iterate
    def main_fuel_flow(self):
        return self.chamber_fuel_flow + self.turbine_mass_flow

    @property
    def turbine_inlet_temperature(self):
        if self.cooling_channel_section.outlet_temperature > self.maximum_wall_temperature:
            warnings.warn('Cooling outlet temperature cannot be higher than maximum wall temperature, cooling flow must be increased manually')
        return self.cooling_channel_section.outlet_temperature


if __name__ == '__main__':
    import arguments as args
    f1 = OpenExpanderCycle(**args.desgin_arguments,**args.base_arguments,**args.oe_arguments)
    print(f1)