from dataclasses import dataclass, field
from functools import cached_property
from math import sqrt, pi
from typing import Optional
from EngineComponents.Abstract.PressureComponent import PressureComponent, NewPressureComponent

#
# @dataclass
# class CombustionChamber(PressureComponent):
#     throat_area: float  # [m2]
#     combustion_chamber_pressure: float  # [Pa]
#     convergent_volume_estimate: float  # [m3]
#     characteristic_length: float  # [m]
#     area_ratio_chamber_throat: float  # [-]
#
#     @cached_property
#     def volume_incl_nozzle_convergent(self):
#         return self.characteristic_length * self.throat_area
#
#     # All subsequent properties of the chamber concern only the cylindrical section, i.e. the actual chamber
#     @cached_property
#     def volume(self):
#         return self.volume_incl_nozzle_convergent - self.convergent_volume_estimate
#
#     @cached_property
#     def geometry_factor(self):
#         return 2
#
#     @cached_property
#     def max_pressure(self):
#         return self.combustion_chamber_pressure
#
#     @cached_property
#     def area(self):
#         return self.area_ratio_chamber_throat * self.throat_area
#
#     @cached_property
#     def radius(self):
#         return sqrt(self.area / pi)
#
#     @cached_property
#     def length(self):
#         return self.volume / self.area
#
#     @cached_property
#     def surface_area(self):
#         return self.radius * 2 * pi * self.length
#

@dataclass
class CombustionChamber(NewPressureComponent):
    throat_area: float  # [m2]
    combustion_chamber_pressure: float  # [Pa]
    convergent_volume_estimate: float  # [m3]
    characteristic_length: float  # [m]
    area_ratio_chamber_throat: float  # [-]
    geometry: str = field(init=False, default='open_cylinder')

    @cached_property
    def volume_incl_nozzle_convergent(self):
        return self.characteristic_length * self.throat_area

    # All subsequent properties of the chamber concern only the cylindrical section, i.e. the actual chamber
    @cached_property
    def volume(self):
        return self.volume_incl_nozzle_convergent - self.convergent_volume_estimate

    @cached_property
    def max_pressure(self):
        return self.combustion_chamber_pressure

    @cached_property
    def area(self):
        return self.area_ratio_chamber_throat * self.throat_area

    @cached_property
    def length(self):
        return self.volume / self.area



