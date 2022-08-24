from dataclasses import dataclass
from math import sqrt, pi
from typing import Optional

from BaseEngineCycle.Structure import Structure
from FlowComponent import FlowComponent

@dataclass
class Injector(Structure, FlowComponent):
    combustion_chamber_pressure: float  # [Pa]
    combustion_chamber_area: float  # [m2]
    propellant_is_gas: Optional[bool] = None
    _pressure_drop_factor = .3  # [-]

    @property
    def inlet_pressure(self):
        return self.inlet_flow_state.pressure

    @property
    def pressure_change(self):
        option1 = self.combustion_chamber_pressure * self._pressure_drop_factor
        # Mota 2008 -> Kesaev and Almeida 2005
        f = .4 if self.propellant_is_gas else .8
        option2 = f * 10E2 * sqrt(10 * self.combustion_chamber_pressure)
        return -option1

    @property
    def thickness(self):
        # Zandbergen -> 4 times pressure drop, ASME code for welded flat caps on internal pressure vessels
        diameter = 2 * sqrt(self.combustion_chamber_area / pi)
        pressure = 4 * self.pressure_drop
        c = 6 * 3.3 / 64
        return diameter * sqrt(c * pressure / self.yield_strength) * self.safety_factor

    @property
    def mass(self):
        return self.combustion_chamber_area * self.thickness * self.material_density