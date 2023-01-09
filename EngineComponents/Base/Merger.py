from dataclasses import dataclass, replace, field
from EngineComponents.Abstract.FlowState import FlowState, DefaultFlowState
from EngineComponents.Abstract.FlowComponent import FlowComponent
import warnings
from math import isclose

from EngineComponents.Base.Splitter import Splitter


@dataclass
class Merger(FlowComponent):
    """ Merges multiple inlet flow states into a single outlet flow state. Sums mass_flows and averages temperatures.
    Pressures and propellants of input flows are expected to be equal.

    :param inlet_flow_states: Multiple flow states that will be combined to create outlet flow state
    """
    inlet_flow_states: tuple[FlowState, ...] = (DefaultFlowState(),)
    is_homogeneous_flows: bool = True
    _warn_pressure: bool = field(init=False, repr=False, default=True)
    inlet_flow_state: FlowState = field(init=False, repr=False, default=DefaultFlowState())

    def __post_init__(self):
        if Merger._warn_pressure:
            self.pressure_check()
        if self.is_homogeneous_flows:
            self.name_check()
        self.set_inlet_flow_state()

    def set_inlet_flow_state(self):
        """Set inlet_flow_state as combination of all inlet states."""
        name = self.inlet_flow_states[0].propellant_name if self.is_homogeneous_flows else 'ChamberGas'
        type = self.inlet_flow_states[0].type if self.is_homogeneous_flows else 'combusted'
        self.inlet_flow_state = replace(self.inlet_flow_states[0],
                                        propellant_name=name,
                                        temperature=self.average_temperature,
                                        mass_flow=self.total_mass_flow,
                                        type=type,)

    def pressure_check(self):
        pressures = [flow_state.pressure for flow_state in self.inlet_flow_states]
        if not all(isclose(pressure, pressures[0]) for pressure in pressures):
            pressures = [f'{flow_state.pressure:.4e}' + " Pa" for flow_state in self.inlet_flow_states]
            warnings.warn(
                'Pressures of incoming streams in Merger are not equal, which could lead to back flow. Ensure '
                'the streams are equal pressure at the inlet of the merger\n'
                'FlowState pressures:\n'
                f'{pressures}')

    def name_check(self):
        names = [flow_state.propellant_name for flow_state in self.inlet_flow_states]
        if not all(name == names[0] for name in names):
            raise ValueError('Inlet streams of Merger must have the same propellant_name')

    @property
    def total_mass_flow(self):
        """"Sums the inlet mass flows"""
        return sum(flow_state.mass_flow for flow_state in self.inlet_flow_states)

    @property
    def average_temperature(self):
        """Calculates the temperature after combining of the inlet flows.

        Assumes the same heat capacity for all flows, thus possible over simplification for non-homogeneous flows.
        """
        return sum(flow_state.mass_flow * flow_state.temperature / self.total_mass_flow
                   for flow_state in self.inlet_flow_states)
