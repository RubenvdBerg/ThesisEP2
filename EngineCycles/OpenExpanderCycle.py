from dataclasses import dataclass, field
from typing import Optional
from EngineCycles.Abstract.OpenCycle import OpenEngineCycle, OpenEngineCycle_DoubleTurbine
from EngineComponents.Base.Merger import Merger
from EngineComponents.Base.Splitter import Splitter
from EngineComponents.Base.Pump import Pump
from EngineCycles.CoolantBleedCycle import BaseCoolantBleedCycle_Mixin


@dataclass
class OpenExpanderCycle_Mixin(BaseCoolantBleedCycle_Mixin):
    """Similar to CoolantBleedCycle, but the fuel flow is split before AND after the cooling channel."""

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
        return self.pre_cooling_splitter.outlet_flow_states['coolant']


@dataclass
class OpenExpanderCycle_DoublePump_Mixin(OpenExpanderCycle_Mixin):
    """Similar to CoolantBleedCycle, but the fuel flow is split before AND after the cooling channel."""

    secondary_fuel_pump_efficiency: Optional[float] = None
    _secondary_fuel_pump_pressure_factor_first_guess: float = .3
    _secondary_fuel_pump_outlet_pressure: float = field(init=False, repr=False, default=None)

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

    @property
    def cooling_inlet_flow_state(self):
        return self.secondary_fuel_pump.outlet_flow_state

    @property
    def fuel_pump_expected_pressure(self):
        return self.combustion_chamber_pressure - self.injector.pressure_change

    @property
    def secondary_fuel_pump_expected_pressure(self):
        return self.fuel_pump_expected_pressure - self.cooling_channel_section.pressure_change

    @property
    def pumps_mass(self):
        return super().pumps_mass + self.secondary_fuel_pump.mass

    @property
    def components_list(self):
        return super().components_list + [
            'Secondary Fuel Pump',
        ]



@dataclass
class OpenExpanderCycle_DoublePump(OpenExpanderCycle_DoublePump_Mixin, OpenEngineCycle):
    @property
    def turbine_mass_flow_initial_guess(self):
        return self.base_mass_flow * .02

    @property
    def fuel_pumps_power_required(self):
        return (self.fuel_pump.power_required + self.secondary_fuel_pump.power_required) / self.shaft_mechanical_efficiency


@dataclass
class OpenExpanderCycle_DoublePumpTurbine(OpenExpanderCycle_DoublePump_Mixin, OpenEngineCycle_DoubleTurbine):

    @property
    def fuel_turbine_mass_flow_initial_guess(self):
        return .01 * self.base_mass_flow

    @property
    def oxidizer_turbine_mass_flow_initial_guess(self):
        return .01 * self.base_mass_flow

    @property
    def fuel_pumps_power_required(self):
        return (self.fuel_pump.power_required + self.secondary_fuel_pump.power_required) / self.shaft_mechanical_efficiency


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


@dataclass
class OpenExpanderCycle_DoubleTurbineSeries(OpenExpanderCycle_DoubleTurbine):
    """GasGeneratorCycle_DoubleTurbine assumes the turbines are in parallel and have separate exhausts. This
    configuration has turbines in series and a single exhaust. The fuel secondary exhaust and associated variables are
    removed."""

    # Variables no longer needed
    fuel_secondary_specific_impulse_quality_factor: float = field(init=False, repr=False)
    fuel_exhaust_expansion_ratio: Optional[float] = field(init=False, repr=False)
    fuel_exhaust_exit_pressure_forced: Optional[float] = field(init=False, repr=False)

    # Properties no longer needed
    @property
    def turbine_splitter(self):
        raise NotImplementedError

    @property
    def fuel_secondary_exhaust(self):
        raise NotImplementedError

    # Adjust iteration
    def turbine_flow_error_larger_than_accuracy(self):
        required = max(self.oxidizer_turbine.mass_flow_required, self.fuel_turbine.mass_flow_required)
        error = abs(required - self.turbine_mass_flow)
        margin = required * self.iteration_accuracy
        return error > margin

    @property
    def turbine_mass_flow(self):
        return max(self._iterative_oxidizer_turbine_mass_flow, self._iterative_fuel_turbine_mass_flow)

    @property
    def verbose_iteration_required(self):
        return max(self.fuel_turbine.mass_flow_required,self.oxidizer_turbine.mass_flow_required)

    # Adjusted inlet flows
    @property
    def fuel_turbine_inlet_flow_state(self):
        return self.post_cooling_splitter.outlet_flow_states['turbine']

    @property
    def oxidizer_turbine_inlet_flow_state(self):
        return self.fuel_turbine.outlet_flow_state

    # Adjusted exhaust thrust
    @property
    def exhaust_total_thrust(self):
        return self.oxidizer_secondary_exhaust.thrust

    @property
    def engine_dry_mass(self):
        return self.feed_system_mass + self.thrust_chamber.mass + self.oxidizer_secondary_exhaust.mass

    @property
    def components_list(self):
        components_list: list = super().components_list
        components_list.remove('Fuel Secondary Exhaust')
        return components_list

