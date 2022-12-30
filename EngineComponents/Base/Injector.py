from dataclasses import dataclass
from math import sqrt, pi
from typing import Optional

from EngineComponents.Abstract.StructuralComponent import StructuralComponent
from EngineComponents.Base.Merger import Merger


@dataclass
class Injector(Merger, StructuralComponent):
    """Injector inherits from Merger to be able to handle both oxidizer en fuel inlet flows"""
    combustion_chamber_pressure: float = 0  # [Pa]
    combustion_chamber_area: float = 0  # [m2]
    pressure_drop_factor: float = 0  # [-]
    is_homogeneous_flows: bool = False

    @property
    def pressure_change(self):
        option1 = self.combustion_chamber_pressure * self.pressure_drop_factor
        # # Mota 2008 -> Kesaev and Almeida 2005
        # f = .4 if self.propellant_is_gas else .8
        # option2 = f * 10E2 * sqrt(10 * self.combustion_chamber_pressure)
        return -option1

    # @property
    # def thickness(self):
    #     # Zandbergen -> 4 times pressure drop, ASME code for welded flat caps on internal pressure vessels
    #     diameter = 2 * sqrt(self.combustion_chamber_area / pi)
    #     pressure = 4 * self.pressure_drop
    #     c = 6 * 3.3 / 64
    #     return diameter * sqrt(c * pressure / self.structure_material.yield_strength) * self.safety_factor

    @property
    def combustion_chamber_radius(self):
        return sqrt(self.combustion_chamber_area / pi)

    @property
    def mass(self):
        return (self.safety_factor * self.structure_material.density * self.combustion_chamber_area
                * sqrt(0.75 * (1 + self.structure_material.poisson_ratio) * self.combustion_chamber_pressure
                       * self.combustion_chamber_radius ** 2 / self.structure_material.yield_strength))
