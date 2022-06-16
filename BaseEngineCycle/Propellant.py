from dataclasses import dataclass
from typing import Literal


@dataclass
class Propellant:
    name: str
    type: Literal['oxidizer', 'fuel']
    mass_flow: float  # [kg/s]
    burn_time: float  # [s]
    density: float  # [kg/m3]
    margin_factor: float  # [-]

    @property
    def mass(self):
        return self.mass_flow * self.burn_time * self.margin_factor

    @property
    def volume(self):
        return self.mass / self.density

    @property
    def volumetric_flow(self):
        return self.mass_flow / self.density



