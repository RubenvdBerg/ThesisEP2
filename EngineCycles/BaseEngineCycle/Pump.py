from dataclasses import dataclass
from typing import Optional
from EngineCycles.BaseEngineCycle.FlowComponent import FlowComponent


@dataclass
class Pump(FlowComponent):
    pressure_increase: float = 0 # [Pa]
    efficiency: float = 0 # [-]
    specific_power: float = 0 # [W/kg]
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
        return self.volumetric_flow_rate * self.pressure_increase / self.efficiency

    @property
    def mass(self):
        return self.power_required / self.specific_power

    @property
    def pressure_change(self):
        return self.pressure_increase
