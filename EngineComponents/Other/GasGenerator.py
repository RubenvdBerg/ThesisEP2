from EngineComponents.Abstract.PressureComponent import PressureComponent, NewPressureComponent
from EngineComponents.Abstract.FlowState import ManualFlowState, DefaultFlowState, FlowState
from dataclasses import dataclass, replace, field
from typing import Optional


@dataclass
class GasGenerator(NewPressureComponent):
    """Handles mass estimation of the Gas Generator and merging of inlet flows.

    The outlet temperature is forced and possibly NOT equal to base_flow_state.temperature!
    mass_mixture_ratio is required for clarity, so it is available as attribute for other components"""
    stay_time: float = 0  # [s]
    base_flow_state: FlowState = DefaultFlowState()
    outlet_temperature: float = 0  # [K]
    oxidizer_inlet_flow_state: FlowState = DefaultFlowState()
    fuel_inlet_flow_state: FlowState = DefaultFlowState()
    mass_mixture_ratio: float = 0  # [-]
    geometry: str = field(init=False, default='sphere')

    @property
    def outlet_mass_flow(self) -> float:
        return sum((self.oxidizer_inlet_flow_state.mass_flow, self.fuel_inlet_flow_state.mass_flow))

    @property
    def outlet_flow_state(self) -> FlowState:
        return replace(self.base_flow_state,
                       mass_flow=self.outlet_mass_flow,
                       temperature=self.outlet_temperature, )

    # Mass Calculation Properties
    @property
    def max_pressure(self):
        return self.base_flow_state.pressure

    @property
    def volume(self):
        return self.stay_time * self.outlet_mass_flow / self.base_flow_state.density
