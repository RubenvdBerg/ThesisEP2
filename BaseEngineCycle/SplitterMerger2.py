from dataclasses import dataclass, replace, field
from typing import Iterable, Optional
from BaseEngineCycle.FlowState import FlowState, DefaultFlowState
from BaseEngineCycle.FlowComponent import BaseFlowComponent, FlowComponent


@dataclass
class Merger(BaseFlowComponent):
    """ Merges multiple inlet flow states into a single outlet flow state. Sums mass_flows and averages temperatures.
    Pressures and propellants of input flows are expected to be equal.

    :param inlet_flow_states: Multiple flow states that will be combined to create outlet flow state
    """
    inlet_flow_states: tuple[FlowState, ...] = (DefaultFlowState(),)

    def __post_init__(self):
        self.pressure_check()
        self.name_check()

    def pressure_check(self):
        pressures = [flow_state.pressure for flow_state in self._inlet_flow_states]
        if not all(pressure == pressures[0] for pressure in pressures):
            flow_states = [flow_state + "\n" for flow_state in self._inlet_flow_states]
            raise ValueError('Pressures of incoming streams are not equal, which could lead to back flow. Ensure the '
                             'streams are equal pressure at the inlet of the merger\n'
                             'FlowStates:\n'
                             f'{flow_states}')

    def name_check(self):
        names = [flow_state.propellant_name for flow_state in self._inlet_flow_states]
        if not all(name == names[0] for name in names):
            raise ValueError('Inlet streams of Merger must have the same propellant_name')

    @property
    def total_mass_flow(self):
        """"Sums the inlet mass flows"""
        return sum(flow_state.mass_flow for flow_state in self._inlet_flow_states)

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


@dataclass
class Splitter(FlowComponent):
    """Splits inlet flow state into multiple outlet flow states.

    :param outlet_mass_flows: If given with N mass flows, creates N+1 outlet flow states with the last one's mass flow
    being equal to the inlet mass flow minus the sum of all N given outlet mass flows
    :param mass_flow_fractions: If given with N mass flow fractions, creates N outlet flow states, each with a
    corresponding fraction of the inlet mass flow
    """

    outlet_mass_flows: Optional[tuple[float, ...]] = None
    mass_flow_fractions: Optional[tuple[float, ...]] = None
    resolved_mass_flows: tuple[float, ...] = field(init=False)
    outlet_flow_names: Optional[tuple[str, ...]] = None

    def __post_init__(self):
        self.split_flows()

    def split_flows(self):
        """Checks values of parameters and delegates handling to the right method"""
        if not ((self.outlet_mass_flows is None) ^ (self.mass_flow_fractions is None)):
            raise ValueError('Exactly one of outlet_mass_flows and mass_flow_fractions must be given')
        elif self.outlet_mass_flows is not None:
            self.resolve_outlet_mass_flows()
        elif self.mass_flow_fractions is not None:
            self.resolve_fractional_outlet_mass_flows()
        self.create_outlet_mass_flows()

    def resolve_outlet_mass_flows(self):
        """Creates N+1 outlet mass flows, the value of the last outlet mass flow taken such that the total sum of outlet
         mass flows equals the inlet mass flow
        """
        if sum(self.outlet_mass_flows) > self.inlet_flow_state.mass_flow:
            raise ValueError('Sum of given mass flows must be less than the inlet mass flow')
        final_mass_flow = self.inlet_mass_flow - sum(self.outlet_mass_flows)
        self.resolved_mass_flows = self.outlet_mass_flows + (final_mass_flow,)

    def resolve_fractional_outlet_mass_flows(self):
        """Creates N outlet mass flows with fractional mass flows of the inlet mass flow"""
        self.resolved_mass_flows = tuple(fraction * self.inlet_mass_flow for fraction in self.mass_flow_fractions)

    def create_outlet_mass_flows(self):
        """"Creates outlet flow states with equal outlet pressures and temperatures (possibly different from inlet state
         using _pressure/_temperature_change"""
        if self.outlet_flow_names is None:
            names = [i for i, _ in enumerate(self.resolved_mass_flows)]
        else:
            if len(self.resolved_mass_flows) != len(self.outlet_flow_names):
                raise ValueError('An equal amount of outlet_flow_names must be given as outlet flow states will be '
                                 'generated, i.e. [len(outlet_mass_flows) + 1] or [len(mass_flow_fractions)] names')
            names = self.outlet_flow_names
        for name, mass_flow in zip(names, self.resolved_mass_flows):
            setattr(self,
                    f'outlet_flow_state_{name}',
                    replace(self.inlet_flow_state,
                            pressure=self.outlet_pressure,
                            temperature=self.outlet_temperature,
                            mass_flow=mass_flow))

    @property
    def outlet_flow_state(self):
        raise ValueError('A Splitter has more than one outlet_flow_state, request outlet_flow_state_[i]/_[name] instead')


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
