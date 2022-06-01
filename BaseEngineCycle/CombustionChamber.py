from dataclasses import dataclass
from functools import cached_property
from math import sqrt, pi
from typing import Optional

from scipy.interpolate import interp1d

from BaseEngineCycle.Structure import Structure


@dataclass
class CombustionChamber(Structure):
    throat_area: float  # [m2]
    combustion_chamber_pressure: float  # [Pa]
    area_ratio_chamber_throat: Optional[float] = None
    propellant_mix: Optional[str] = None
    characteristic_length: Optional[float] = None

    def __post_init__(self):
        if self.characteristic_length is None:
            try:
                self.characteristic_length = self.characteristic_length_options[self.propellant_mix]
            except KeyError:
                raise KeyError(
                    f'The specified propellant mix[{self.propellant_mix}] does not have a default characteristic length'
                    f' for the combustion chamber, specify one manually or select a propellant mix that does ['
                    f'{self.characteristic_length_options.keys}]')
        if self.area_ratio_chamber_throat is None:
            # Humble 1995 p.222
            throat_diameter = 2 * sqrt(self.throat_area / pi)
            self.area_ratio_chamber_throat = (8.0 * throat_diameter ** 2.4 + 1.25)

    @cached_property
    def characteristic_length_options(self):
        return {'LOX/GH2': 0.635, 'LOX/LH2': 0.89, 'LOX/RP1': 1.145}

    @property
    def volume(self):
        return self.characteristic_length * self.throat_area

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
    def mass(self):
        return self.material_density * self.safety_factor / self.yield_strength * (2 * self.volume
                                                                                   * self.combustion_chamber_pressure)


@dataclass
class ComplexCombustionChamber(CombustionChamber):
    max_wall_temperature: Optional[float] = None  # [K]
    material_name: Optional[str] = None

    def __post_init__(self):
        if self.yield_strength is None:
            data = self.yield_strength_data[self.material_name]
            self.yield_strength = self.yield_strength_from_temp(
                temp_list=data['temps'], sigma_list=data['sigmas'], max_temperature=self.max_wall_temperature
            )

    @staticmethod
    def yield_strength_from_temp(temp_list: list, sigma_list: list, max_temperature: float) -> float:
        def yield_strength_function(temp: float) -> float:
            f = interp1d(temp_list, sigma_list)
            return f(temp) * 1e6

        return yield_strength_function(max_temperature)

    @cached_property
    def yield_strength_data(self):
        # Temps in [K], sigmas [MPa]
        return {'inconel600':
                {'temps': [
                    297,
                    600,
                    800,
                    900,
                    1000,
                    1050,
                    1100,
                    1150,
                    1200,
                    1300,
                    1373],
                    'sigmas': [
                        733,
                        689,
                        617,
                        468,
                        273,
                        212,
                        154,
                        113,
                        79,
                        50,
                        27]}}


@dataclass
class Injector(Structure):
    combustion_chamber_pressure: float  # [Pa]
    combustion_chamber_area: float  # [m2]
    propellant_is_gas: bool

    @property
    def pressure_drop(self):
        # Mota 2008 -> Kesaev and Almeida 2005
        f = .4 if self.propellant_is_gas else .8
        return f * 10E2 * sqrt(10 * self.combustion_chamber_pressure)

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