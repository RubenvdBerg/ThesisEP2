from dataclasses import dataclass
from typing import Optional

from EngineComponents.Abstract.FlowComponent import FlowComponent


@dataclass
class BatteryCooler(FlowComponent):
    """Component that adjusts its outlet flow to be equal to the coolant flow required. Used to iterate until flows
    through pump and battery are matching in the EP Cycle"""
    power_heat_loss: float = 0  # [W]
    outlet_pressure_required: float = 0  # [Pa]
    coolant_allowable_temperature_change: float = 0  # [K]
    coolant_specific_heat_capacity: Optional[float] = None  # [J/(kgK)]

    def __post_init__(self):
        if self.coolant_specific_heat_capacity is None:
            self.coolant_specific_heat_capacity = self.inlet_flow_state.specific_heat_capacity

    @property
    def mass_flow_change(self):
        return self.coolant_flow_required - self.inlet_mass_flow

    @property
    def pressure_change(self):
        return self.outlet_pressure_required - self.inlet_pressure

    @property
    def temperature_change(self):
        return self.coolant_allowable_temperature_change

    @property
    def coolant_flow_required(self):
        return self.power_heat_loss / (self.coolant_specific_heat_capacity * self.coolant_allowable_temperature_change)