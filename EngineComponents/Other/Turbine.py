from dataclasses import dataclass
from typing import Optional
from EngineComponents.Abstract.FlowComponent import FlowComponent


@dataclass
class Turbine(FlowComponent):
    power_required: float = 0  # [W]
    efficiency: float = 0  # [-]
    pressure_ratio: Optional[float] = None  # [-]
    outlet_pressure_forced: Optional[float] = None  # [Pa]

    def __post_init__(self):
        self.resolve_pressure_ratio_choice()

    def resolve_pressure_ratio_choice(self):
        if not ((self.pressure_ratio is None) ^ (self.outlet_pressure_forced is None)):
            raise ValueError('Both or neither the pressure_ratio and the outlet_pressure of the turbine are provided. '
                             'Provide one and only one')
        elif self.pressure_ratio is None:
            self.pressure_ratio = self.inlet_pressure / self.outlet_pressure_forced

    @property
    def mass_flow_required(self):
        cp = self.inlet_flow_state.specific_heat_capacity
        y = self.inlet_flow_state.heat_capacity_ratio
        return (self.power_required / (self.efficiency * cp * self.inlet_temperature
                                       * (1 - self.pressure_ratio ** ((1 - y) / y))))

    @property
    def temperature_change(self):
        if self.mass_flow == 0:
            return 0
        return -1 * self.power_required / (self.mass_flow * self.inlet_flow_state.specific_heat_capacity)

    @property
    def pressure_change(self):
        return self.inlet_pressure / self.pressure_ratio - self.inlet_pressure
