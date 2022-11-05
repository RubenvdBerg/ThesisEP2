import warnings
from dataclasses import dataclass
from math import pi
from typing import Optional

import numpy as np

from EngineCycles.BaseEngineCycle.Structure import PressureStructure
from EngineCycles.BaseEngineCycle.FlowComponent import FlowComponent


@dataclass
class Tank(FlowComponent, PressureStructure):
    max_acceleration: float = 0 # [m/s2]
    ullage_factor: float = 0  # [-]
    propellant_volume: float = 0  # [m3]
    pressurant_tank_volume: Optional[float] = None  # [m3]
    propellant_density: Optional[float] = None

    @property
    def _propellant_density(self):
        if self.propellant_density is None:
            return self.inlet_flow_state.density
        else:
            warnings.warn(f'Manual {self.inlet_flow_state.type} density provided:[{self.propellant_density:.2f}], estimated density:[{self.inlet_flow_state.density:.2f}]')
            return self.propellant_density

    @property
    def initial_pressure(self):
        return self.inlet_flow_state.pressure

    @property
    def radius(self):
        return (3 * self.volume / (4 * pi)) ** (1 / 3)

    @property
    def unused_volume(self):
        ullage_volume = (self.ullage_factor - 1) * self.propellant_volume
        if self.pressurant_tank_volume is None:
            return ullage_volume
        else:
            return self.pressurant_tank_volume + ullage_volume

    @property
    def cap_height(self):
        r= self.radius
        a = pi / 3
        b = -pi * r
        c = 0
        d = self.unused_volume
        sols = np.roots(np.array([a,b,c,d]))
        sols = [sol for sol in sols if r > sol > 0]
        return min(sols)

    @property
    def initial_head(self):
        return 2 * self.radius - self.cap_height

    @property
    def total_lower_pressure(self):
        return self.initial_pressure + self._propellant_density * self.max_acceleration * self.initial_head

    @property
    def total_upper_pressure(self):
        return self.initial_pressure + self._propellant_density * self.max_acceleration * (self.initial_head -
                                                                                                self.radius)

    @property
    def maximum_expected_operating_pressure(self):
        return (self.total_upper_pressure + self.total_lower_pressure) / 2

    @property
    def volume(self):
        if self.pressurant_tank_volume is not None:
            return self.propellant_volume * self.ullage_factor + self.pressurant_tank_volume
        else:
            warnings.warn(
                f'No pressurant tank volume was given. Assumed pressurant tank is not submerged in {self.inlet_flow_state.propellant_name} tank')
            return self.propellant_volume * self.ullage_factor