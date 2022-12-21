import warnings
from dataclasses import dataclass, field, replace
from functools import cached_property
from typing import Optional

from EngineCycles.BaseOpenCycle.OpenCycle import OpenEngineCycle, OpenEngineCycle_DoubleTurbine
from EngineCycles.BaseEngineCycle.Merger import Merger
from EngineCycles.BaseEngineCycle.Splitter import Splitter
from EngineCycles.BaseEngineCycle.FlowState import FlowState
from EngineCycles.BaseEngineCycle.Cooling import CoolingChannelSection
from EngineCycles.BaseEngineCycle.Pump import Pump


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


@dataclass
class OpenExpanderCycle_Mixin(BaseCoolantBleedCycle_Mixin):
    """Similar to Coolant Bleed Cycle, but the fuel flow is split before and after the cooling channel.

     Just like the CoolantBleedCycle the fuel flow is split into a chamber flow and turbine flow, but less pressure is
     required and the turbine flow will be hotter. Which is also the risk of this flow, the mass flow required by the
     turbine might be less than the mass flow required for cooling, consequently the mass flow to the cooling and
     turbine needs to be increased, which is bad because this flow is dumped overboard at low specific impulse after
     the turbine, so less flow split off is better"""

    secondary_fuel_pump_efficiency: Optional[float] = None
    _secondary_fuel_pump_pressure_factor_first_guess: float = .3
    _secondary_fuel_pump_outlet_pressure: float = field(init=False, repr=False, default=None)

    @property
    def required_coolant_mass_flow(self):
        """Set the required coolant mass flow to max(turbine_flow, min_req_cool_flow). Additionally, if engine is still
        initializing (m_flow_cool = 0) or if the actual coolant flow will be equal to the minimum required flow; the
        temp calc is ignored, which leads to setting the cooling_channel_section.outlet_temperature to
        maximum_coolant_outlet_temperature.
        """
        if self.minimum_required_coolant_mass_flow == 0 or self.minimum_required_coolant_mass_flow > self.turbine_mass_flow:
            self._is_temp_calc_needed = False
            return self.minimum_required_coolant_mass_flow
        else:
            self._is_temp_calc_needed = True
            return self.turbine_mass_flow

    def calc_pump_outlet_pressures(self):
        super().calc_pump_outlet_pressures()
        self._secondary_fuel_pump_outlet_pressure = self.secondary_fuel_pump_expected_pressure

    @property
    def secondary_fuel_pump_outlet_pressure(self):
        # Use estimate that is based on chamber pressure, until calc_pump_outlet_pressures has run
        if self._secondary_fuel_pump_outlet_pressure:
            return self._secondary_fuel_pump_outlet_pressure
        else:
            return self.combustion_chamber_pressure * self._secondary_fuel_pump_pressure_factor_first_guess + self.fuel_pump_outlet_pressure

    @property
    def secondary_fuel_pump(self):
        eta2 = self.secondary_fuel_pump_efficiency
        eta = self.fuel_pump_efficiency if eta2 is None else eta2
        return Pump(inlet_flow_state=self.pre_cooling_splitter.outlet_flow_state_coolant,
                    expected_outlet_pressure=self.secondary_fuel_pump_outlet_pressure,
                    efficiency=eta,
                    specific_power=self.fuel_pump_specific_power)

    # New components required to split the flow before AND after the cooling channels
    @property
    def pre_cooling_splitter(self):
        """Split fuel flow into required coolant flow and "rest" flow, which should be equal to primary chamber flow."""
        return Splitter(inlet_flow_state=self.fuel_pump.outlet_flow_state,
                        required_outlet_mass_flows=(self.required_coolant_mass_flow,),
                        outlet_flow_names=('coolant', 'chamber'))

    @property
    def post_cooling_splitter(self):
        """Split coolant (fuel) flow into second chamber flow and "rest" flow, which should be equal to turbine flow."""
        chamber_flow1 = self.pre_cooling_splitter.outlet_flow_state_chamber.mass_flow
        return Splitter(inlet_flow_state=self.cooling_channel_section.outlet_flow_state,
                        required_outlet_mass_flows=(self.chamber_fuel_flow - chamber_flow1,),
                        outlet_flow_names=('chamber', 'turbine'))

    @property
    def pre_injection_merger(self):
        """Merge primary and secondary chamber fuel flows"""
        return Merger(inlet_flow_states=(self.post_cooling_splitter.outlet_flow_state_chamber,
                                         self.pre_cooling_splitter.outlet_flow_state_chamber))

    @property
    def turbine_inlet_flow_state(self):
        return self.post_cooling_splitter.outlet_flow_states['turbine']

    @property
    def injector_inlet_flow_states(self):
        return self.pre_injection_merger.outlet_flow_state, self.oxidizer_pump.outlet_flow_state

    @property
    def cooling_inlet_flow_state(self):
        return self.secondary_fuel_pump.outlet_flow_state

    @property
    def pumps_power_required(self):
        return (
                           self.fuel_pump.power_required + self.secondary_fuel_pump.power_required + self.oxidizer_pump.power_required) / self.shaft_mechanical_efficiency

    @property
    def fuel_pump_expected_pressure(self):
        return self.combustion_chamber_pressure - self.injector.pressure_change

    @property
    def secondary_fuel_pump_expected_pressure(self):
        return self.fuel_pump_expected_pressure - self.cooling_channel_section.pressure_change


@dataclass
class OpenExpanderCycle(OpenExpanderCycle_Mixin, OpenEngineCycle):
    @property
    def turbine_mass_flow_initial_guess(self):
        return self.base_mass_flow * .02


@dataclass
class OpenExpanderCycle_DoubleTurbine(OpenExpanderCycle_Mixin, OpenEngineCycle_DoubleTurbine):

    @property
    def fuel_turbine_mass_flow_initial_guess(self):
        return .01 * self.base_mass_flow

    @property
    def oxidizer_turbine_mass_flow_initial_guess(self):
        return .01 * self.base_mass_flow


if __name__ == '__main__':
    import arguments as args

    f1 = CoolantBleedCycle(**args.desgin_arguments, **args.base_arguments, **args.oe_arguments)
    print(f1.cooling_channel_section.throat_wall_temperature(f1.maximum_wall_temperature))
    print(f1)
