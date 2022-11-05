from dataclasses import dataclass, replace
from EngineCycles.BaseEngineCycle.FlowState import FlowState, DefaultFlowState
from EngineCycles.BaseEngineCycle.FlowComponent import BaseFlowComponent
import warnings

from EngineCycles.BaseEngineCycle.Splitter import Splitter


@dataclass
class Merger(BaseFlowComponent):
    """ Merges multiple inlet flow states into a single outlet flow state. Sums mass_flows and averages temperatures.
    Pressures and propellants of input flows are expected to be equal.

    :param inlet_flow_states: Multiple flow states that will be combined to create outlet flow state
    """
    inlet_flow_states: tuple[FlowState, ...] = (DefaultFlowState(),)
    is_homogeneous_flows: bool = True

    def __post_init__(self):
        self.pressure_check()
        if self.is_homogeneous_flows:
            self.name_check()

    def pressure_check(self):
        pressures = [flow_state.pressure for flow_state in self.inlet_flow_states]
        if not all(pressure == pressures[0] for pressure in pressures):
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
        """Calculates the temperature after combining of the inlet flows"""
        return sum(flow_state.mass_flow * flow_state.temperature / self.total_mass_flow
                   for flow_state in self.inlet_flow_states)

    @property
    def _inlet_flow_state(self) -> FlowState:
        """Overwrites Parents property to also account for average temperature"""
        return replace(self.inlet_flow_states[0],
                       temperature=self.average_temperature,
                       mass_flow=self.total_mass_flow)


if __name__ == '__main__':
    kwargs = {'propellant_name': 'RP-1', 'pressure': 1E6, 'type': 'fuel'}
    s1 = FlowState(temperature=300, mass_flow=10, **kwargs)
    s2 = FlowState(temperature=400, mass_flow=20, **kwargs)
    merger = Merger(inlet_flow_states=(s1, s2))
    print(merger.outlet_flow_state)

    m1 = Splitter(inlet_flow_state=s1, outlet_mass_flows=(8.,))
    m2 = Splitter(inlet_flow_state=s1, mass_flow_fractions=(.2, .8))
    print(m2.outlet_flow_state_0)
    print(m2.outlet_flow_state_1)
    print(m1.outlet_flow_state_0)
    print(m1.outlet_flow_state_1)
