from dataclasses import dataclass
from typing import Optional
from EngineComponents.Abstract.FlowComponent import FlowComponent


@dataclass
class Pump(FlowComponent):
    expected_outlet_pressure: float = 0  # [Pa]
    # pressure_increase: float = 0  # [Pa]
    efficiency: float = 0  # [-]
    specific_power: float = 0  # [W/kg]
    propellant_density: Optional[float] = None

    @property
    def _propellant_density(self):
        if self.propellant_density is None:
            return self.inlet_flow_state.density
        else:
            return self.propellant_density

    @property
    def volumetric_flow_rate(self):
        return self.inlet_mass_flow / self._propellant_density

    @property
    def power_required(self):
        if self.pressure_change < 0:
            raise ValueError('Negative pressure change required over pump, increase combustion pressure or decrease tank pressure')
        return self.volumetric_flow_rate * self.pressure_change / self.efficiency

    @property
    def mass(self):
        return self.power_required / self.specific_power

    @property
    def pressure_change(self):
        return self.expected_outlet_pressure - self.inlet_pressure
