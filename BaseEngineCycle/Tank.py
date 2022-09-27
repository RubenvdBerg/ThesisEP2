import warnings
from dataclasses import dataclass
from math import pi
from typing import Optional
from BaseEngineCycle.Propellant import Propellant
from BaseEngineCycle.Structure import PressureStructure
from BaseEngineCycle.FlowComponent import FlowComponent


@dataclass
class Tank(FlowComponent, PressureStructure):
    max_acceleration: float = 0 # [m/s2]
    ullage_factor: float = 0  # [-]
    propellant_volume: float = 0  # [m3]
    pressurant_tank_volume: Optional[float] = None  # [m3]

    @property
    def initial_pressure(self):
        return self.inlet_flow_state.pressure

    @property
    def radius(self):
        return (3 * self.volume / (4 * pi)) ** (1 / 3)

    @property
    def initial_head(self):
        unused_volume = (self.ullage_factor - 1) * self.propellant_volume
        rest_height = (unused_volume * 3 / (2 * pi)) ** (1 / 3)
        return 2 * self.radius - rest_height

    @property
    def total_lower_pressure(self):
        return self.initial_pressure + self.inlet_flow_state.density * self.max_acceleration * self.initial_head

    @property
    def total_upper_pressure(self):
        return self.initial_pressure + self.inlet_flow_state.density * self.max_acceleration * (self.initial_head -
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

# @dataclass
# class Tank(Structure):
#     propellant: Propellant
#     max_acceleration: float  # [m/s2]
#     ullage_factor: float  # [-]
#     initial_pressure: float  # [Pa]
#     kwak_fix_cycle_type: str
#     pressurant_tank_volume: Optional[float] = None  # [m3]
#     kwak_fix: bool = False
#
#     @property
#     def radius(self):
#         return (3 * self.volume / (4 * pi)) ** (1 / 3)
#
#     @property
#     def initial_head(self):
#         unused_volume = (self.ullage_factor - 1) * self.propellant_volume
#         rest_height = (unused_volume * 3 / (2 * pi)) ** (1 / 3)
#         if self.kwak_fix:
#             if self.kwak_fix_cycle_type == 'ep':
#                 if self.propellant.type == 'oxidizer':
#                     return 1.91839449096392
#                 elif self.propellant.type == 'fuel':
#                     return 1.65560478870526
#             elif self.kwak_fix_cycle_type == 'gg':
#                 if self.propellant.type == 'oxidizer':
#                     return 1.92539045861846
#                 elif self.propellant.type == 'fuel':
#                     return 1.71033897378923
#         return 2 * self.radius - rest_height
#
#     @property
#     def total_lower_pressure(self):
#         return self.initial_pressure + self.propellant.density * self.max_acceleration * self.initial_head
#
#     @property
#     def total_upper_pressure(self):
#         return self.initial_pressure + self.propellant.density * self.max_acceleration * (self.initial_head -
#                                                                                           self.radius)
#
#     @property
#     def mass(self):  # Spherical Tanks
#         return (self.safety_factor * 3 * self.material_density / (4 * self.yield_strength) *
#                 self.volume * (self.total_upper_pressure + self.total_lower_pressure))
#
#     @property
#     def volume(self):
#         if self.pressurant_tank_volume is not None:
#             return self.propellant_volume * self.ullage_factor + self.pressurant_tank_volume
#         else:
#             warnings.warn(
#                 f'No pressurant tank volume was given. Assumed pressurant tank is not submerged in {self.propellant.type} tank')
#             return self.propellant_volume * self.ullage_factor
