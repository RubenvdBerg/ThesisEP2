from dataclasses import dataclass
from EngineComponents.Abstract.StructuralComponent import StructuralComponent


@dataclass
class PressureComponent(StructuralComponent):

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
        return (self.geometry_factor * self.safety_factor * self.structure_material.density
                * self.maximum_expected_operating_pressure * self.volume / self.structure_material.yield_strength)


