from dataclasses import dataclass
from functools import cached_property
from math import pi, sqrt, cos, sin, tan, asin
from typing import Optional
import scipy.optimize
from numpy import array


def get_chamber_throat_area_ratio_estimate(throat_area: float):
    # Humble 1995 p.222
    throat_diameter_in_cm = 2 * sqrt(throat_area / pi) * 1e2
    area_ratio_chamber_throat = (8.0 * throat_diameter_in_cm ** -.6 + 1.25)
    return area_ratio_chamber_throat

@dataclass
class Nozzle:
    # conv, div = convergent and divergent sections of the nozzle respectively
    throat_area: float  # [m2]
    expansion_ratio: float  # [-]
    conv_chamber_bend_ratio: float  # [-]
    conv_throat_bend_ratio: float  # [-]
    conv_half_angle: float  # [rad]
    div_throat_half_angle: float  # [rad]
    div_longi_throat_radius: Optional[float] = None  # [rad]
    area_ratio_chamber_throat: Optional[float] = None

    def __post_init__(self):
        if self.conv_half_angle > pi / 2 or self.conv_half_angle < 0:
            raise ValueError(
                f'Half angle of the convergent must be between 0 and \u03C0/2 (Given:{self.conv_half_angle:.3f})')
        if self.area_ratio_chamber_throat is None:
            self.area_ratio_chamber_throat = get_chamber_throat_area_ratio_estimate(self.throat_area)
        # [MODERN ENGINEERING FOR DESIGN OF LIQUID-PROPELLANT ROCKET ENGINES, Huzel&Huang 1992, p.76, fig. 4-15]
        # radius of the nozzle after the throat curve and distance between throat and end of throat curve
        self.div_radius_p = self.throat_radius + (1 - cos(self.div_throat_half_angle)) * self.div_longi_throat_radius
        self.div_length_p = self.div_longi_throat_radius * sin(self.div_throat_half_angle)

    @property
    def exit_area(self):
        return self.throat_area * self.expansion_ratio

    @property
    def chamber_area(self):
        return self.area_ratio_chamber_throat * self.throat_area

    @cached_property
    def throat_radius(self):
        return sqrt(self.throat_area / pi)

    @property
    def exit_radius(self):
        return sqrt(self.exit_area / pi)

    @property
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
                             f'please change the bend ratios or divergence half angle')
        return self.conv_throat_long_radius * sin(self.conv_half_angle) + z / tan(self.conv_half_angle)

    @cached_property
    def conv_length(self):
        return self.conv_length_p + self.conv_chamber_long_radius * sin(self.conv_half_angle)

    @property
    def conv_volume_estimate(self):
        r1 = self.throat_radius
        r2 = self.chamber_radius
        h = self.conv_length
        return pi / 3 * h * (r1**2 + r1*r2 + r2**2)

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
        return float(radius)

    def div_radius(self, distance_from_throat: float) -> float:
        # Distance from throat positive towards nozzle exit
        if distance_from_throat > self.div_length or distance_from_throat < 0:
            raise ValueError(
                f'Nozzle diameter cannot be calculated before the throat [<0m] or after the nozzle exit [>{self.div_length:.4E}m].')
        if distance_from_throat < self.div_length_p:
            alpha = asin(distance_from_throat / self.div_longi_throat_radius) / 2
            div_radius = self.throat_radius + distance_from_throat * tan(alpha)
        else:
            div_radius = self.div_radius_after_throat_curve(distance_from_throat)
        return float(div_radius)

    def div_radius_after_throat_curve(self, distance_from_throat:float) -> float:
        raise NotImplementedError('Abstract Class Nozzle does not have a defined divergent radius after the throat curve')

    def get_radius(self, distance_from_throat: float) -> float:
        if distance_from_throat < 0:
            return self.conv_radius(-distance_from_throat)
        elif distance_from_throat > 0:
            return self.div_radius(distance_from_throat)
        else:
            return self.throat_radius


@dataclass
class BellNozzle(Nozzle):
    # conv, div are the convergent and divergent sections of the nozzle respectively
    div_exit_half_angle: float = 0  # [rad]

    def __post_init__(self):
        if self.div_longi_throat_radius is None:
            # As suggested by Huzel&Huang
            self.div_longi_throat_radius = .382 * self.throat_radius
        super().__post_init__()
        # [MODERN ENGINEERING FOR DESIGN OF LIQUID-PROPELLANT ROCKET ENGINES, Huzel&Huang 1992, p.76, fig. 4-15]
        # parabolic equation parameters
        tan_th = tan(pi/2 - self.div_throat_half_angle)
        tan_ex = tan(pi/2 - self.div_exit_half_angle)
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
    def div_length(self):
        a = self.div_a
        b = self.div_b
        c = self.div_c
        y = self.exit_radius
        return a * y**2 + b * y + c

    def div_radius_after_throat_curve(self, distance_from_throat:float) -> float:
        def func(x):
            a = float(self.div_a * x ** 2 + self.div_b * x + self.div_c - distance_from_throat)
            return array([a], dtype=float)

        x0 = float(self.throat_radius + (self.exit_radius - self.throat_radius) * (
                distance_from_throat / self.div_length))
        div_radius, *_ = scipy.optimize.fsolve(func, array([x0], dtype=float))
        return float(div_radius)


    # def div_radius(self, distance_from_throat: float) -> float:
    #     # Distance from throat positive towards nozzle exit
    #     if distance_from_throat > self.div_length or distance_from_throat < 0:
    #         raise ValueError(
    #             f'Nozzle diameter cannot be calculated before the throat [<0m] or after the nozzle exit [>{self.div_length:.4E}m].')
    #     if distance_from_throat < self.div_length_p:
    #         div_radius = self.div_radius_throat_curve(distance_from_throat)
    #     else:
    #         def func(x):
    #             a = float(self.div_a * x ** 2 + self.div_b * x + self.div_c - distance_from_throat)
    #             return array([a], dtype=float)
    #
    #         x0 = float(self.throat_radius + (self.exit_radius - self.throat_radius) * (
    #                     distance_from_throat / self.div_length))
    #         div_radius = scipy.optimize.fsolve(func, array([x0], dtype=float))
    #     return float(div_radius)


@dataclass
class ConicalNozzle(Nozzle):

    def __post_init__(self):
        if self.div_longi_throat_radius is None:
            self.div_longi_throat_radius = .5 * self.throat_radius
        super().__post_init__()

    @property
    def div_length(self):
        e = self.expansion_ratio
        r_t = self.throat_radius
        r_u = self.div_longi_throat_radius
        theta = self.div_throat_half_angle
        part = (e**.5 - 1) * r_t + r_u * (1/cos(theta) - 1)
        return part / tan(theta)

    def div_radius_after_throat_curve(self, distance_from_throat: float) -> float:
        return self.div_radius_p + (distance_from_throat - self.div_length_p) * tan(self.div_throat_half_angle)


if __name__ == '__main__':
    from math import radians
    nozzle = BellNozzle(throat_area=(.1**2)*pi,expansion_ratio=50, conv_throat_bend_ratio=.8, conv_half_angle=radians(20), conv_chamber_bend_ratio=1., div_throat_half_angle=radians(30), div_exit_half_angle=radians(5))
    print()