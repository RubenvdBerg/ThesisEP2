import warnings
from dataclasses import dataclass

from BaseOpenCycle.OpenCycle import OpenEngineCycle
from BaseEngineCycle.SplitterMerger2 import Splitter

@dataclass
class CoolantBleedCycle(OpenEngineCycle):

    @property
    def turbine_flow_splitter(self):
        """Splits flow into required chamber flow and "rest flow" which should be the turbine flow"""
        return Splitter(inlet_flow_state=self.cooling_channel_section.outlet_flow_state,
                        required_outlet_mass_flows=(self.chamber_fuel_flow,),
                        outlet_flow_names=('chamber', 'turbine'))


    @property
    def get_injector_inlet_flow_state_fuel(self):
        return self.turbine_flow_splitter.outlet_flow_state_chamber

    @property
    def default_turbine_flow_check_state(self):
        """Flow provided here should be the outlet flow state of the component before the turbine and will be compared
        to the iteratively found turbine inlet flow state and will warn you if these are not found to be equal"""
        return self.turbine_flow_splitter.outlet_flow_state_turbine

    @property
    def turbine_mass_flow_initial_guess(self):
        return self.base_mass_flow * .02

    @property
    def main_fuel_flow(self):
        return self.chamber_fuel_flow + self.turbine_mass_flow

    @property
    def turbine_inlet_temperature(self):
        if self.cooling_channel_section.outlet_temperature > self.maximum_wall_temperature:
            warnings.warn('Cooling outlet temperature cannot be higher than maximum wall temperature, cooling flow must be increased manually')
        return self.cooling_channel_section.outlet_temperature

@dataclass
class OpenExpanderCycle(CoolantBleedCycle):
    """Similar to Coolant Bleed Cycle, but the fuel flow is split before the cooling channel.

     Just like the CoolantBleedCycle the fuel flow is split into a chamber flow and turbine flow, but less pressure is
     required and the turbine flow will be hotter. Which is also the risk of this flow, the mass flow required by the
     turbine might be less than the mass flow required for cooling, consequently the mass flow to the cooling and
     turbine needs to be increased, which is bad because this flow is dumped overboard at low specific impulse after
     the turbine, so less flow split off is better"""

    @property
    def turbine_flow_splitter(self):
        """In this configuration the turbine splitter is BEFORE the cooling section"""
        return Splitter(inlet_flow_state=self.fuel_pump.outlet_flow_state,
                        required_outlet_mass_flows=(self.chamber_fuel_flow,),
                        outlet_flow_names=('chamber', 'coolant'))


    @property
    def cooling_inlet_flow_state(self):
        return self.turbine_flow_splitter.outlet_flow_state_coolant

    @property
    def default_turbine_flow_check_state(self):
        return self.cooling_channel_section.outlet_flow_state




if __name__ == '__main__':
    import arguments as args
    f1 = CoolantBleedCycle(**args.desgin_arguments, **args.base_arguments, **args.oe_arguments)
    print(f1.cooling_channel_section.throat_wall_temperature(f1.maximum_wall_temperature))
    print(f1)