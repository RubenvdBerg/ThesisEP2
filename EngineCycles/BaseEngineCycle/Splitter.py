from dataclasses import dataclass, field, replace
from typing import Optional

from EngineCycles.BaseEngineCycle.FlowComponent import FlowComponent


@dataclass
class Splitter(FlowComponent):
    """Splits inlet flow state into multiple outlet flow states.

    :param outlet_mass_flows: If given with N mass flows, creates N+1 outlet flow states with the last one's mass flow
    being equal to the inlet mass flow minus the sum of all N given outlet mass flows
    :param mass_flow_fractions: If given with N mass flow fractions, creates N outlet flow states, each with a
    corresponding fraction of the inlet mass flow
    """

    required_outlet_mass_flows: Optional[tuple[float, ...]] = None
    mass_flow_fractions: Optional[tuple[float, ...]] = None
    outlet_flow_names: Optional[tuple[str, ...]] = None
    resolved_mass_flows: tuple[float, ...] = field(init=False)
    outlet_flow_states: dict = field(init=False, default_factory=dict)

    def __post_init__(self):
        self.split_flows()

    def split_flows(self):
        """Checks values of parameters and delegates handling to the right method"""
        if not ((self.required_outlet_mass_flows is None) ^ (self.mass_flow_fractions is None)):
            raise ValueError('Exactly one of outlet_mass_flows and mass_flow_fractions must be given')
        elif self.required_outlet_mass_flows is not None:
            self.resolve_outlet_mass_flows()
        elif self.mass_flow_fractions is not None:
            self.resolve_fractional_outlet_mass_flows()
        self.create_outlet_flow_states()

    def resolve_outlet_mass_flows(self):
        """Creates N+1 outlet mass flows, the value of the last outlet mass flow taken such that the total sum of outlet
         mass flows equals the inlet mass flow
        """
        if sum(self.required_outlet_mass_flows) > self.inlet_flow_state.mass_flow:
            raise ValueError('Sum of given mass flows must be less than the inlet mass flow')
        final_mass_flow = self.inlet_mass_flow - sum(self.required_outlet_mass_flows)
        self.resolved_mass_flows = self.required_outlet_mass_flows + (final_mass_flow,)

    def resolve_fractional_outlet_mass_flows(self):
        """Creates N outlet mass flows with fractional mass flows of the inlet mass flow"""
        self.resolved_mass_flows = tuple(fraction * self.inlet_mass_flow for fraction in self.mass_flow_fractions)

    def create_outlet_flow_states(self):
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
            state = replace(self.inlet_flow_state,
                            pressure=self.outlet_pressure,
                            temperature=self.outlet_temperature,
                            mass_flow=mass_flow)
            self.outlet_flow_states[name] = state
            setattr(self,
                    f'outlet_flow_state_{name}',
                    state)

    @property
    def outlet_flow_state(self):
        raise ValueError(
            'A Splitter has more than one outlet_flow_state, request outlet_flow_state_{name} or outlet_flow_states["name"] instead')