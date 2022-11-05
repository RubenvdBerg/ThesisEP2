import warnings
from dataclasses import dataclass, field, replace
from functools import cached_property
from typing import Optional

from EngineCycles.BaseOpenCycle.OpenCycle import OpenEngineCycle
from EngineCycles.BaseEngineCycle.Merger import Merger
from EngineCycles.BaseEngineCycle.Splitter import Splitter
from EngineCycles.BaseEngineCycle.FlowState import FlowState
from EngineCycles.BaseEngineCycle.Cooling import CoolingChannelSection
from EngineCycles.BaseEngineCycle.Pump import Pump


@dataclass
class CoolantBleedCycle(OpenEngineCycle):
    @property
    def post_cooling_splitter(self):
        """Splits flow into required chamber flow and "rest flow" which should be the turbine flow"""
        return Splitter(inlet_flow_state=self.cooling_channel_section.outlet_flow_state,
                        required_outlet_mass_flows=(self.chamber_fuel_flow,),
                        outlet_flow_names=('chamber', 'turbine'))

    @property
    def injector_inlet_flow_state_fuel(self):
        return self.post_cooling_splitter.outlet_flow_state_chamber

    @property
    def default_turbine_flow_check_state(self):
        """Flow provided here should be the outlet flow state of the component before the turbine and will be compared
        to the iteratively found turbine inlet flow state and will warn you if these are not found to be equal"""
        return self.post_cooling_splitter.outlet_flow_state_turbine

    @property
    def turbine_mass_flow_initial_guess(self):
        return self.base_mass_flow * .02

    @property
    def main_fuel_flow(self):
        return self.chamber_fuel_flow + self.turbine_mass_flow

    @property
    def cooling_channel_section(self):
        """Overwrite parent cooling to use maximum turbine temp instead of maximum wall temp"""
        ccs = super().cooling_channel_section
        ccs.maximum_outlet_temperature = self.turbine_maximum_temperature
        return ccs

    @property
    def turbine_inlet_temperature(self):
        if self.cooling_channel_section.outlet_temperature > self.maximum_wall_temperature:
            warnings.warn(
                'Cooling outlet temperature cannot be higher than maximum wall temperature, cooling flow must be increased manually')
        return self.cooling_channel_section.outlet_temperature


@dataclass
class OpenExpanderCycle(CoolantBleedCycle):
    """Similar to Coolant Bleed Cycle, but the fuel flow is split before the cooling channel.

     Just like the CoolantBleedCycle the fuel flow is split into a chamber flow and turbine flow, but less pressure is
     required and the turbine flow will be hotter. Which is also the risk of this flow, the mass flow required by the
     turbine might be less than the mass flow required for cooling, consequently the mass flow to the cooling and
     turbine needs to be increased, which is bad because this flow is dumped overboard at low specific impulse after
     the turbine, so less flow split off is better"""

    secondary_fuel_pump_pressure_change_factor: float = .2
    secondary_fuel_pump_efficiency: Optional[float] = None
    _is_temp_calc_needed: bool = field(init=False, default=True)

    @property
    def coolant_base_state(self):
        """This state is used instead of cooling_inlet_flow_state to prevent recursion
        (mass flow is dependent on splitters and splitters need this input).
        """
        return FlowState(propellant_name=self.fuel_initial_flow_state.propellant_name,
                         temperature=self.fuel_initial_temperature,
                         pressure=self.fuel_initial_pressure + self.delta_p_fuel_pump + self.delta_p_secondary_fuel_pump,
                         mass_flow=0,
                         type='fuel')

    @cached_property
    def minimum_required_coolant_mass_flow(self):
        """Determine minimum coolant flow required to keep outlet temp below turbine maximum."""
        # Copy of CCS with new state
        ccs = CoolingChannelSection(inlet_flow_state=self.coolant_base_state,
                                    total_heat_transfer=self.heat_transfer_section.total_heat_transfer,
                                    maximum_outlet_temperature=self.turbine_maximum_temperature,
                                    combustion_chamber_pressure=self.combustion_chamber_pressure,
                                    pressure_drop=self.cooling_section_pressure_drop,
                                    verbose=self.verbose, )
        return ccs.min_mass_flow

    @property
    def required_coolant_mass_flow(self):
        if self.minimum_required_coolant_mass_flow < self.turbine_mass_flow:
            self._is_temp_calc_needed = True
            return self.turbine_mass_flow
        else:
            self._is_temp_calc_needed = False
            return self.minimum_required_coolant_mass_flow

    @property
    def cooling_outlet_flow_state(self):
        if self._is_temp_calc_needed:
            return self.cooling_channel_section.outlet_flow_state
        else:
            # Dummy CoolingChannelSection used to find the outlet pressure
            dummy_ccs = CoolingChannelSection(inlet_flow_state=self.cooling_inlet_flow_state,
                                              total_heat_transfer=0,
                                              maximum_outlet_temperature=self.turbine_maximum_temperature,
                                              pressure_drop=self.cooling_section_pressure_drop,
                                              combustion_chamber_pressure=self.combustion_chamber_pressure,
                                              verbose=self.verbose)

            return replace(self.cooling_inlet_flow_state,
                           temperature=self.maximum_wall_temperature,
                           pressure=dummy_ccs.outlet_pressure)

    # @property
    # def cooling_outlet_flow_state(self):
    #     return self.cooling_channel_section.outlet_flow_state

    # New components required to split the flow before AND after the cooling channels
    @property
    def pre_cooling_splitter(self):
        """Split fuel flow into required coolant flow and "rest" flow, which should be equal to primary chamber flow."""
        return Splitter(inlet_flow_state=self.fuel_pump.outlet_flow_state,
                        required_outlet_mass_flows=(self.required_coolant_mass_flow,),
                        outlet_flow_names=('coolant', 'chamber'))

    @property
    def delta_p_secondary_fuel_pump(self):
        return self.secondary_fuel_pump_pressure_change_factor * self.combustion_chamber_pressure

    @property
    def secondary_fuel_pump(self):
        eta2 = self.secondary_fuel_pump_efficiency
        eta = self.fuel_pump_efficiency if eta2 is None else eta2
        return Pump(inlet_flow_state=self.pre_cooling_splitter.outlet_flow_state_coolant,
                    pressure_increase=self.delta_p_secondary_fuel_pump,
                    efficiency=eta,
                    specific_power=self.fuel_pump_specific_power)

    @property
    def post_cooling_splitter(self):
        """Split coolant (fuel) flow into second chamber flow and "rest" flow, which should be equal to turbine flow."""
        chamber_flow1 = self.pre_cooling_splitter.outlet_flow_state_chamber.mass_flow
        return Splitter(inlet_flow_state=self.cooling_outlet_flow_state,
                        required_outlet_mass_flows=(self.chamber_fuel_flow - chamber_flow1,),
                        outlet_flow_names=('chamber', 'turbine'))

    @property
    def pre_injection_merger(self):
        """Merge primary and secondary chamber fuel flows"""
        return Merger(inlet_flow_states=(self.post_cooling_splitter.outlet_flow_state_chamber,
                                         self.pre_cooling_splitter.outlet_flow_state_chamber))

    @property
    def injector_inlet_flow_states(self):
        return self.pre_injection_merger.outlet_flow_state, self.oxidizer_pump.outlet_flow_state

    @property
    def cooling_inlet_flow_state(self):
        return self.secondary_fuel_pump.outlet_flow_state


if __name__ == '__main__':
    import arguments as args

    f1 = CoolantBleedCycle(**args.desgin_arguments, **args.base_arguments, **args.oe_arguments)
    print(f1.cooling_channel_section.throat_wall_temperature(f1.maximum_wall_temperature))
    print(f1)
