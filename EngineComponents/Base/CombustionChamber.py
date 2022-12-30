from dataclasses import dataclass
from functools import cached_property
from math import sqrt, pi
from typing import Optional
from EngineComponents.Abstract.PressureComponent import PressureComponent
from EngineFunctions.EmpiricalRelations import get_chamber_throat_area_ratio_estimate


@dataclass
class CombustionChamber(PressureComponent):
    throat_area: float  # [m2]
    combustion_chamber_pressure: float  # [Pa]
    convergent_volume_estimate: float  # [m3]
    characteristic_length: float
    area_ratio_chamber_throat: Optional[float] = None
    propellant_mix: Optional[str] = None

    def __post_init__(self):
        if self.area_ratio_chamber_throat is None:
            self.area_ratio_chamber_throat = get_chamber_throat_area_ratio_estimate(self.throat_area)

    @property
    def volume_incl_nozzle_convergent(self):
        return self.characteristic_length * self.throat_area

    # All subsequent properties of the chamber concern only the cylindrical section, i.e. the actual chamber
    @property
    def volume(self):
        return self.volume_incl_nozzle_convergent - self.convergent_volume_estimate

    @property
    def area(self):
        return self.area_ratio_chamber_throat * self.throat_area

    @property
    def radius(self):
        return sqrt(self.area / pi)

    @property
    def length(self):
        return self.volume / self.area

    @property
    def geometry_factor(self):
        return 2

    @property
    def maximum_expected_operating_pressure(self):
        return self.combustion_chamber_pressure

    @property
    def surface_area(self):
        return self.radius * 2 * pi * self.length


