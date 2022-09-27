from dataclasses import dataclass

@dataclass
class Structure:
    material_density: float  # [kg/m3]
    safety_factor: float  # [-]
    yield_strength: float  # [Pa]


@dataclass
class PressureStructure(Structure):
    @property
    def volume(self):
        raise NotImplementedError

    @property
    def maximum_expected_operating_pressure(self):
        raise NotImplementedError

    @property
    def geometry_factor(self):
        """Geometrical factor for thin walled pressure vessels,
        r = radius, l = length, v = volume
        3/2                     for spheres,
        2 + (2r/l)              for flat capped cylinders
        3 - (pi * l * r**2 / v) for cylinders with hemispherical caps"""
        return 3/2

    @property
    def mass(self):
        return (self.geometry_factor * self.safety_factor * self.material_density
                * self.maximum_expected_operating_pressure * self.volume)
