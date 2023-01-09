from dataclasses import dataclass
from functools import cached_property
from math import pi, sqrt, cos, sin, tan, asin
from numpy import array, linspace, interp, concatenate, trapz
from EngineComponents.Abstract.Material import Material
from EngineComponents.Abstract.StructuralComponent import StructuralComponent


@dataclass
class Nozzle(StructuralComponent):
    # ConicalNozzle
    # conv, div = convergent and divergent sections of the nozzle respectively
    throat_area: float  # [m2]
    expansion_ratio: float  # [-]
    conv_chamber_bend_ratio: float  # [-]
    conv_throat_bend_ratio: float  # [-]
    conv_half_angle: float  # [rad]
    div_throat_half_angle: float  # [rad]
    chamber_pressure: float  # [Pa]
    area_ratio_chamber_throat: float  # [-]

    _interpol_num: float = 40

    def __post_init__(self):
        if self.conv_half_angle > pi / 2 or self.conv_half_angle < 0:
            raise ValueError(
                f'Half angle of the convergent must be between 0 and \u03C0/2 (Given:{self.conv_half_angle:.3f})')
        # [MODERN ENGINEERING FOR DESIGN OF LIQUID-PROPELLANT ROCKET ENGINES, Huzel&Huang 1992, p.76, fig. 4-15]
        # radius of the nozzle after the throat curve and distance between throat and end of throat curve
        self.div_radius_p = self.throat_radius + (1 - cos(self.div_throat_half_angle)) * self.div_longi_throat_radius
        self.div_length_p = self.div_longi_throat_radius * sin(self.div_throat_half_angle)

        self.set_radius_func()

    @cached_property
    def end_interpol_distances(self):
        return array([self.div_length])

    def set_radius_func(self):
        x_cc_bend = linspace(-self.conv_length, -self.conv_length_p, self._interpol_num)
        length_q = self.conv_throat_long_radius * sin(self.conv_half_angle)
        x_throat = linspace(-length_q, self.div_length_p, self._interpol_num)
        self.x_distances = concatenate((x_cc_bend, x_throat, self.end_interpol_distances))
        self.y_radii = array([self.get_radius_original(distance) for distance in self.x_distances])
        # self.surface_area = trapz(y_radii * 2 * pi, x=x_distances)

    def get_radius(self, distance_from_throat):
        return interp(distance_from_throat, self.x_distances, self.y_radii)

    @cached_property
    def exit_area(self):
        return self.throat_area * self.expansion_ratio

    @cached_property
    def chamber_area(self):
        return self.area_ratio_chamber_throat * self.throat_area

    @cached_property
    def throat_radius(self):
        return sqrt(self.throat_area / pi)

    @cached_property
    def exit_radius(self):
        return sqrt(self.exit_area / pi)

    @cached_property
    def chamber_radius(self):
        return sqrt(self.chamber_area / pi)

    @cached_property
    def conv_throat_long_radius(self):
        # Convergent longitudinal radius at throat [m]
        return self.conv_throat_bend_ratio * self.throat_radius

    @cached_property
    def conv_chamber_long_radius(self):
        # Convergent longitudinal radius at connection to combustion chamber
        return self.conv_chamber_bend_ratio * self.chamber_radius

    @cached_property
    def conv_length_p(self):
        z = (self.chamber_radius - self.throat_radius
             - (self.conv_throat_long_radius + self.conv_chamber_long_radius) * (1 - cos(self.conv_half_angle)))
        if z < 0:
            raise ValueError(f'The length of the straight part of the convergent is calculated to be less than zero, '
                             f'please change the bend ratios or chamber/throat area ratio.')
        return self.conv_throat_long_radius * sin(self.conv_half_angle) + z / tan(self.conv_half_angle)

    @cached_property
    def conv_length(self):
        return self.conv_length_p + self.conv_chamber_long_radius * sin(self.conv_half_angle)

    @cached_property
    def conv_volume_estimate(self):
        r1 = self.throat_radius
        r2 = self.chamber_radius
        h = self.conv_length
        return pi / 3 * h * (r1 ** 2 + r1 * r2 + r2 ** 2)

    @cached_property
    def total_length(self):
        return self.conv_length + self.div_length

    def get_conv_radius(self, distance_from_throat: float) -> float:
        # Zandbergen 2017 Lecture Notes p.66-69
        # Distance from throat positive towards chamber
        r_u = self.conv_throat_long_radius
        r_a = self.conv_chamber_long_radius
        length_q = r_u * sin(self.conv_half_angle)
        radius_q = self.throat_radius + r_u * (1 - cos(self.conv_half_angle))
        radius_p = self.chamber_radius - r_a * (1 - cos(self.conv_half_angle))
        length_p = self.conv_length_p
        if distance_from_throat > self.conv_length or distance_from_throat < 0:
            raise ValueError(
                f'Radius of the convergent cannot be calculated above its length ({self.conv_length:.4e} or downstream of the throat')
        elif distance_from_throat > length_p:
            distance = self.conv_length - distance_from_throat
            alpha = asin(distance / r_a) / 2
            radius = self.chamber_radius - distance * tan(alpha)
        elif distance_from_throat > length_q:
            radius = radius_q + (distance_from_throat - length_q) * tan(self.conv_half_angle)
        else:
            alpha = asin(distance_from_throat / r_u) / 2
            radius = self.throat_radius + distance_from_throat * tan(alpha)
        return float(radius)

    @cached_property
    def div_longi_throat_radius(self):
        # Divergent longitudinal radius at throat [m]
        return .5 * self.throat_radius

    @cached_property
    def div_length(self):
        e = self.expansion_ratio
        r_t = self.throat_radius
        r_u = self.div_longi_throat_radius
        theta = self.div_throat_half_angle
        part = (e ** .5 - 1) * r_t + r_u * (1 / cos(theta) - 1)
        return part / tan(theta)

    def get_div_radius_after_throat_curve(self, distance_from_throat: float) -> float:
        return self.div_radius_p + (distance_from_throat - self.div_length_p) * tan(self.div_throat_half_angle)

    def get_div_radius(self, distance_from_throat: float) -> float:
        # Distance from throat positive towards nozzle exit
        if distance_from_throat > self.div_length or distance_from_throat < 0:
            raise ValueError(
                f'Radius of the nozzle divergent cannot be calculated before the throat [<0m] or after the nozzle exit [>{self.div_length:.4E}m].')
        if distance_from_throat < self.div_length_p:
            alpha = asin(distance_from_throat / self.div_longi_throat_radius) / 2
            div_radius = self.throat_radius + distance_from_throat * tan(alpha)
        else:
            div_radius = self.get_div_radius_after_throat_curve(distance_from_throat)
        return float(div_radius)

    def get_radius_original(self, distance_from_throat: float) -> float:
        if distance_from_throat < 0:
            return self.get_conv_radius(-distance_from_throat)
        elif distance_from_throat > 0:
            return self.get_div_radius(distance_from_throat)
        else:
            return self.throat_radius

    # Mass estimation properties
    @cached_property
    def conv_surface_area(self):
        # Simplified as curved surface of a frustrum
        return pi * self.conv_length * (self.chamber_radius + self.throat_radius)

    @cached_property
    def div_surface_area(self):
        return pi * self.div_length * (self.exit_radius + self.throat_radius)

    @cached_property
    def surface_area(self):
        return self.conv_surface_area + self.div_surface_area

    @cached_property
    def wall_thickness(self):
        return self.chamber_pressure * self.chamber_radius / self.structure_material.yield_strength

    @cached_property
    def mass(self):
        return self.safety_factor * self.surface_area * self.wall_thickness * self.structure_material.density


@dataclass
class BellNozzle(Nozzle):
    # conv, div are the convergent and divergent sections of the nozzle respectively
    div_exit_half_angle: float = 0  # [rad]

    def __post_init__(self):
        super().__post_init__()
        # [MODERN ENGINEERING FOR DESIGN OF LIQUID-PROPELLANT ROCKET ENGINES, Huzel&Huang 1992, p.76, fig. 4-15]
        # parabolic equation parameters
        tan_th = tan(pi / 2 - self.div_throat_half_angle)
        tan_ex = tan(pi / 2 - self.div_exit_half_angle)
        self.div_a = ((tan_ex - tan_th) / (2 * (self.exit_radius - self.div_radius_p)))
        self.div_b = tan_th - 2 * self.div_a * self.div_radius_p
        self.div_c = self.div_length_p - self.div_a * self.div_radius_p ** 2 - self.div_b * self.div_radius_p

        cot_th = 1 / tan(self.div_throat_half_angle)
        cot_ex = 1 / tan(self.div_exit_half_angle)
        self.a = (cot_ex - cot_th) / (2 * (self.exit_radius - self.div_radius_p))
        self.b = cot_ex - 2 * self.a * self.exit_radius
        self.b2 = cot_th - 2 * self.a * self.div_radius_p
        self.c = self.div_length_p - self.a * (self.div_radius_p ** 2) - self.b * self.div_radius_p
        print()

    @property
    def div_longi_throat_radius(self):
        # As suggested by Huzel&Huang 1992
        return .382 * self.throat_radius

    @property
    def div_length(self):
        a = self.div_a
        b = self.div_b
        c = self.div_c
        y = self.exit_radius
        return a * y ** 2 + b * y + c

    def get_div_radius_after_throat_curve(self, distance_from_throat: float) -> float:
        def func(x):
            a = float(self.div_a * x ** 2 + self.div_b * x + self.div_c - distance_from_throat)
            return array([a], dtype=float)

        x0 = float(self.throat_radius + (self.exit_radius - self.throat_radius) * (
                distance_from_throat / self.div_length))
        div_radius, *_ = scipy.optimize.fsolve(func, array([x0], dtype=float))
        return float(div_radius)

    @property
    def end_interpol_distances(self):
        return linspace(self.div_length_p, self.div_length, self._interpol_num)
