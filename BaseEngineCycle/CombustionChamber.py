from dataclasses import dataclass
from functools import cached_property
from math import sqrt, pi
from typing import Optional

from scipy.interpolate import interp1d

from BaseEngineCycle.Structure import PressureStructure
from BaseEngineCycle.Nozzle import get_chamber_throat_area_ratio_estimate
from BaseEngineCycle.FlowComponent import FlowComponent

@dataclass
class CombustionChamber(PressureStructure):
    throat_area: float  # [m2]
    combustion_chamber_pressure: float  # [Pa]
    convergent_volume_estimate: float  # [m3]
    area_ratio_chamber_throat: Optional[float] = None
    propellant_mix: Optional[str] = None
    characteristic_length: Optional[float] = None
    verbose: bool = True

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
            self.area_ratio_chamber_throat = get_chamber_throat_area_ratio_estimate(self.throat_area)
        if self.convergent_volume_estimate is None:
            self.convergent_volume_estimate = 0
            if self.verbose:
                warnings.warn(
                    'No estimate for the volume of the convergent given. Calculating combustion chamber length'
                    ' based solely on the volume of the cylindrical section without convergent nozzle section.')

    @cached_property
    def characteristic_length_options(self):
        return {'LOX/GH2': 0.635, 'LOX/LH2': 0.89, 'LOX/RP1': 1.145, 'LOX/LCH4': 1.45}

    @property
    def volume(self):  # Of the cylindrical section
        return self.characteristic_length * self.throat_area - self.convergent_volume_estimate

    @property
    def volume_incl_convergent(self):
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
    def geometry_factor(self):
        return 2 * (1 + self.radius / self.length)

    @property
    def maximum_expected_operating_pressure(self):
        return self.combustion_chamber_pressure

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
        super().__post_init__()

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