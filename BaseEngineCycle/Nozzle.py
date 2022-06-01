from dataclasses import dataclass
from functools import cached_property
from math import pi, sqrt, cos, sin, tan, asin

import scipy.optimize
from numpy import array


@dataclass
class Nozzle:
    # conv, div = convergent and divergent sections of the nozzle respectively
    throat_area: float  # [m2]
    area_ratio: float  # [-]
    chamber_radius: float  # [m]
    conv_chamber_bend_ratio: float  # [-]
    conv_throat_bend_ratio: float  # [-]
    conv_half_angle: float  # [rad]

    def __post_init__(self):
        if self.conv_half_angle > pi / 2 or self.conv_half_angle < 0:
            raise ValueError(
                f'Half angle of the convergent must be between 0 and \u03C0/2 (Given:{self.conv_half_angle:.3f})')

    @cached_property
    def exit_area(self):
        return self.throat_area * self.area_ratio

    @cached_property
    def throat_radius(self):
        return sqrt(self.throat_area / pi)

    @cached_property
    def exit_radius(self):
        return sqrt(self.exit_area / pi)

    @cached_property
    def conv_throat_long_radius(self):
        # Convergent longitudinal radius at throat [m]
        return self.conv_throat_bend_ratio * self.throat_radius

    @cached_property
    def conv_chamber_long_radius(self):
        # Convergent longitudinal radius at connection to combustion chamber
        return self.conv_chamber_bend_ratio * self.chamber_radius

    @property
    def conv_length_p(self):
        z = (self.chamber_radius - self.throat_radius
             - (self.conv_throat_long_radius + self.conv_chamber_long_radius) * (1 - cos(self.conv_half_angle)))
        if z < 0:
            raise ValueError(f'The length of the straight part of the convergent is calculated to be less than zero, '
                             f'please change the bend ratios or divergence half angle')
        return self.conv_throat_long_radius * sin(self.conv_half_angle) + z / tan(self.conv_half_angle)

    @property
    def conv_length(self):
        return self.conv_length_p + self.conv_chamber_long_radius * sin(self.conv_half_angle)

    @property
    def div_length(self):
        raise NotImplementedError('Abstract base class Nozzle does not have a divergent length')

    @property
    def total_length(self):
        return self.conv_length + self.div_length

    def conv_radius(self, distance_from_throat: float) -> float:
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
        return radius

    def div_radius(self, distance_from_throat: float) -> float:
        # Distance from throat positive towards nozzle exit
        raise NotImplementedError('Abstract base class Nozzle does not have a divergent radius')

    def get_radius(self, distance_from_throat: float) -> float:
        if distance_from_throat < 0:
            return self.conv_radius(-distance_from_throat)
        elif distance_from_throat > 0:
            return self.div_radius(distance_from_throat)
        else:
            return self.throat_radius


@dataclass
class BellNozzle(Nozzle):
    # conv, div = convergent and divergent sections of the nozzle respectively
    div_throat_half_angle: float  # [rad]
    div_exit_half_angle: float  # [rad]

    def __post_init__(self):
        super().__post_init__()
        # [MODERN ENGINEERING FOR DESIGN OF LIQUID-PROPELLANT ROCKET ENGINES, Huzel&Huang 1992, p.76, fig. 4-15]
        # radius of the nozzle after the throat curve and distance between throat and end of throat curve
        self.div_radius_p = 1.382 * self.throat_radius - 0.382 * self.throat_radius * cos(self.div_throat_half_angle)
        self.div_length_p = 0.382 * self.throat_radius * sin(self.div_throat_half_angle)
        # parabolic equation parameters
        self.div_a = tan(pi / 2 - self.div_exit_half_angle) - tan(pi / 2 - self.div_throat_half_angle) / (
                2 * (self.exit_radius - self.div_radius_p))
        self.div_b = tan(pi / 2 - self.div_throat_half_angle) - 2 * self.div_a * self.div_radius_p
        self.div_c = self.div_length_p - self.div_a * self.div_radius_p ** 2 - self.div_b * self.div_radius_p

    @property
    def div_length(self):
        return self.div_a * self.exit_radius ** 2 + self.div_b * self.exit_radius + self.div_c

    def div_radius(self, distance_from_throat: float) -> float:
        # Distance from throat positive towards nozzle exit
        if distance_from_throat > self.div_length or distance_from_throat < 0:
            raise ValueError(
                f'Nozzle diameter cannot be calculated before the throat [<0m] or after the nozzle exit [>{self.div_length:.4E}m].')
        if distance_from_throat < self.div_length_p:
            alpha = asin(distance_from_throat / (0.382 * self.throat_radius)) / 2
            div_radius = self.throat_radius + distance_from_throat * tan(alpha)
        else:
            def func(x):
                return array(float(self.div_a * x ** 2 + self.div_b * x + self.div_c - distance_from_throat))

            x0 = array(float(
                self.throat_radius + (self.exit_radius - self.throat_radius) * (
                        distance_from_throat / self.div_length)))
            div_radius = scipy.optimize.fsolve(func, x0)
        return div_radius


@dataclass
class ConicalNozzle(Nozzle):
    # conv, div = convergent and divergent sections of the nozzle respectively
    divergent_half_angle: float  # [rad]

    @property
    def div_length(self):
        raise NotImplementedError

    def div_radius(self, distance_from_throat: float) -> float:
        raise NotImplementedError