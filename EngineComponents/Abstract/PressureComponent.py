import warnings
from typing import Literal
from dataclasses import dataclass, field
from EngineComponents.Abstract.StructuralComponent import StructuralComponent
from math import pi

@dataclass
class PressureComponent(StructuralComponent):
    @property
    def volume(self):
        raise NotImplementedError

    @property
    def max_pressure(self):
        raise NotImplementedError

    @property
    def geometry_factor(self):
        """Geometrical factor for thin walled pressure vessels,
        r = radius, l = length, v = volume
        3/2                     for spheres,
        2                       for open-ended cylinders,
        2 + (2r/l)              for flat-capped cylinders,
        3 - (pi * l * r**2 / v) for cylinders with hemispherical caps"""
        return 3/2

    @property
    def mass(self):
        return (self.geometry_factor * self.safety_factor * self.structure_material.density
                * self.max_pressure * self.volume / self.structure_material.yield_strength)


@dataclass
class NewPressureComponent(StructuralComponent):
    geometry: Literal['sphere', 'open_cylinder', 'closed_cylinder']

    @property
    def max_pressure(self):
        raise NotImplementedError

    @property
    def volume(self):
        raise NotImplementedError

    @property
    def length(self):  # Only for cylinders
        raise NotImplementedError

    @property
    def radius(self):
        if self.geometry == 'sphere':
            return (3 * self.volume / 4 / pi)**(1 / 3)
        else:
            return (self.volume / (pi * self.length)) ** 0.5

    @property
    def surface_area(self):
        if self.geometry == 'sphere':
            return 4 * pi * self.radius**2
        elif self.geometry == 'open_cylinder':
            return 2 * pi * self.radius * self.length
        elif self.geometry == 'closed_cylinder':
            return 2 * pi * self.radius * (self.length + self.radius)

    @property
    def geometry_factor(self):
        return 2 if self.geometry == 'cylinder' else 1

    @property
    def thickness(self):
        calculated_thickness = self.geometry_factor * self.safety_factor * self.max_pressure * self.radius / self.structure_material.yield_strength
        if calculated_thickness < self.structure_material.minimal_thickness:
            warnings.warn('Calculated thickness of pressure component is smaller than minimum thickness of material. The latter is used')
            return self.structure_material.minimal_thickness
        return calculated_thickness

    @property
    def mass(self):
        return self.surface_area * self.thickness * self.structure_material.density
