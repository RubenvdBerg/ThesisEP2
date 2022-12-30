from dataclasses import dataclass
from functools import cached_property
from typing import Optional

from scipy.interpolate import interp1d

from EngineComponents.Base.CombustionChamber import CombustionChamber


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