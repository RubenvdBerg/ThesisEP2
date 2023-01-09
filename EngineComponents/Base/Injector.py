from dataclasses import dataclass, field
from math import sqrt, pi
from typing import Optional

from EngineComponents.Abstract.StructuralComponent import StructuralComponent
from EngineComponents.Base.Merger import Merger


@dataclass
class Injector(Merger, StructuralComponent):
    """Injector inherits from Merger to be able to handle both oxidizer en fuel inlet flows"""
    combustion_chamber_pressure: float = 0  # [Pa]
    combustion_chamber_area: float = 0  # [m2]
    pressure_drop: float = 0  # [-]
    is_homogeneous_flows: bool = field(init=False, repr=False, default=False)

    @property
    def pressure_change(self):
        return -self.pressure_drop

    @property
    def combustion_chamber_radius(self):
        return sqrt(self.combustion_chamber_area / pi)

    @property
    def mass(self):
        return (self.safety_factor * self.structure_material.density * self.combustion_chamber_area
                * sqrt(0.75 * (1 + self.structure_material.poisson_ratio) * self.combustion_chamber_pressure
                       * self.combustion_chamber_radius ** 2 / self.structure_material.yield_strength))
