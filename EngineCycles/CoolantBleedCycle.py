from dataclasses import dataclass
from EngineComponents.Base.Splitter import Splitter
from EngineCycles.Abstract.OpenCycle import OpenEngineCycle


@dataclass
class BaseCoolantBleedCycle_Mixin:
    @property
    def main_fuel_flow(self):
        return self.chamber_fuel_flow + self.turbine_mass_flow

    @property
    def maximum_coolant_outlet_temperature(self):
        return min(self.turbine_maximum_temperature, self.maximum_wall_temperature)


@dataclass
class CoolantBleedCycle_Mixin(BaseCoolantBleedCycle_Mixin):

    @property
    def post_cooling_splitter(self):
        """Splits flow into required chamber flow and "rest flow" which should be the turbine flow"""
        return Splitter(inlet_flow_state=self.cooling_channel_section.outlet_flow_state,
                        required_outlet_mass_flows=(self.chamber_fuel_flow,),
                        outlet_flow_names=('chamber', 'turbine'))

    @property
    def turbine_inlet_flow_state(self):
        return self.post_cooling_splitter.outlet_flow_states['turbine']

    @property
    def injector_inlet_flow_states(self):
        return self.post_cooling_splitter.outlet_flow_states['chamber'], self.oxidizer_pump.outlet_flow_state


@dataclass
class CoolantBleedCycle2_Mixin(BaseCoolantBleedCycle_Mixin):
    """Same as CoolantBleedCycle but flow is split before coolingsection"""

    @property
    def pre_cooling_splitter(self):
        """Split fuel flow into required coolant flow and "rest" flow, which should be equal to primary chamber flow."""
        return Splitter(inlet_flow_state=self.fuel_pump.outlet_flow_state,
                        required_outlet_mass_flows=(self.required_coolant_mass_flow,),
                        outlet_flow_names=('cooling', 'chamber'))

    @property
    def cooling_inlet_flow_state(self):
        return self.pre_cooling_splitter.outlet_flow_states['cooling']

    @property
    def injector_inlet_flow_states(self):
        return self.pre_cooling_splitter.outlet_flow_states['chamber'], self.oxidizer_pump.outlet_flow_state

    @property
    def turbine_inlet_flow_state(self):
        return self.cooling_channel_section.outlet_flow_state

    @property
    def fuel_pump_expected_pressure(self):
        dp1 = self.combustion_chamber_pressure - self.injector.pressure_change
        dp2 = self.combustion_chamber_pressure - self.cooling_channel_section.pressure_change
        return max(dp1, dp2)


@dataclass
class CoolantBleedCycle(CoolantBleedCycle_Mixin, OpenEngineCycle):

    @property
    def turbine_mass_flow_initial_guess(self):
        return self.base_mass_flow * .02