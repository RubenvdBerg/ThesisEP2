from dataclasses import dataclass


@dataclass
class Structure:
    material_density: float  # [kg/m3]
    safety_factor: float  # [-]
    yield_strength: float  # [Pa]