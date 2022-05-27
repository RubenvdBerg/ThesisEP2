import warnings
from math import pi, radians, cos, sin, tan, asin, sqrt
from typing import Callable, Optional

import scipy.integrate
import scipy.optimize
from matplotlib import pyplot as plt
from numpy import array, isclose, linspace
from dataclasses import dataclass

from scipy.interpolate import interp1d

from BaseEngineCycle.base_engine_cycle import Structure
from irt import get_area_ratio

@dataclass
class ThrustChamber:
    nozzle: Nozzle
    combustion_chamber: CombustionChamber


class ThrustChamber:
    def __init__(self, nozzle: type(Nozzle), combustion_chamber: type(CombustionChamber), heat_capacity_ratio: float):
        self.nozzle = nozzle
        self.cc = combustion_chamber
        self.heat_capacity_ratio = heat_capacity_ratio

    def get_radius(self, distance_from_throat):
        if distance_from_throat < -self.cc.total_length:
            raise ValueError(
                f'Radius of thrust chamber cannot be calculated before the injector plate ({self.cc.total_length:.4e} m before the throat)')
        elif distance_from_throat < -self.cc.convergent_length:
            return self.cc.radius
        elif distance_from_throat < 0:
            return self.cc.convergent_radius(-distance_from_throat)
        else:
            return self.nozzle.radius(distance_from_throat)

    def get_mach(self, distance_from_throat):
        radius = self.get_radius(distance_from_throat)
        area_ratio = pi * radius ** 2 / self.throat_area
        x0 = .1 if distance_from_throat < 0 else 10

        def func(mach_number):
            return array(float(get_area_ratio(mach_number, self.heat_capacity_ratio) - area_ratio))
        mach = scipy.optimize.fsolve(func, array(float(x0)))
        check = isclose(func(mach), [0.0])
        if not check:
            warnings.warn('Implicitly calculated Mach number is not within tolerances')
        return float(mach)

    def show_mach(self, **kwargs):
        self.distance_plot(self.get_mach, 'Mach [-]', **kwargs)

    def show_contour(self, **kwargs):
        self.distance_plot(self.get_radius, 'Radius [m]', **kwargs)

    @property
    def min_distance_from_throat(self):
        return -self.cc.total_length

    @property
    def max_distance_from_throat(self):
        return self.nozzle.length

    @property
    def throat_distance_tuple(self):
        return -self.cc.total_length, self.nozzle.length

    @property
    def surface(self):
        result = scipy.integrate.quad(lambda x: self.get_radius(x) * 2 * pi, *self.throat_distance_tuple)
        return result[0]

    @property
    def throat_area(self):
        return self.nozzle.a_t

    def distance_plot(self, func: Callable, ylabel: str, num=300, ytick_function: Optional[Callable] = None):
        distances = linspace(float(-self.cc.total_length), self.nozzle.length, num)
        values = [func(distance) for distance in distances]
        fig, ax = plt.subplots()
        ax.plot(distances, values)
        ax.set_ylabel(ylabel)
        ax.set_xlabel('Distance from throat [m]')
        if ytick_function is not None:
            ticks = ax.get_yticks().tolist()
            ax.set_yticks(ticks)
            ax.set_yticklabels([ytick_function(x) for x in ticks])
        plt.show()


class Nozzle:
    def __init__(self, throat_area, nozzle_type, area_ratio, throat_half_angle, exit_half_angle,
                 angles_in_rad=True):
        self.a_t = throat_area
        self.eps = area_ratio
        self.type = nozzle_type
        self.half_angle_th = throat_half_angle if angles_in_rad else radians(throat_half_angle)
        self.half_angle_ex = exit_half_angle if angles_in_rad else radians(exit_half_angle)
        if self.type == 'bell':
            # [MODERN ENGINEERING FOR DESIGN OF LIQUID-PROPELLANT ROCKET ENGINES, Huzel&Huang 1992, p.76, fig. 4-15]
            # radius of the nozzle after the throat curve and distance between throat and end of throat curve
            self.radius_p = 1.382 * self.throat_radius - 0.382 * self.throat_radius * cos(self.half_angle_th)
            self.length_p = 0.382 * self.throat_radius * sin(self.half_angle_th)
            # parabolic equation parameters
            self.a = tan(pi / 2 - self.half_angle_ex) - tan(pi / 2 - self.half_angle_th) / (
                    2 * (self.exit_radius - self.radius_p))
            self.b = tan(pi / 2 - self.half_angle_th) - 2 * self.a * self.radius_p
            self.c = self.length_p - self.a * self.radius_p ** 2 - self.b * self.radius_p
        elif self.type == 'conical':
            raise NotImplementedError
        else:
            raise NotImplementedError('Only bell nozzles have been implemented at this stage')

    @property
    def length(self):
        if self.type == 'bell':
            return self.a * self.exit_radius ** 2 + self.b * self.exit_radius + self.c
        elif self.type == 'conical':
            raise NotImplementedError

    def radius(self, distance_from_throat: float):
        if distance_from_throat > self.length or distance_from_throat < 0:
            raise ValueError(
                f'Nozzle diameter cannot be calculated before the throat [<0m] or after the nozzle exit [>{self.length:.4E}m].')
        if distance_from_throat < self.length_p:
            alpha = asin(distance_from_throat / (0.382 * self.throat_radius)) / 2
            radius = self.throat_radius + distance_from_throat * tan(alpha)
        else:
            def func(x):
                return array(float(self.a*x**2+self.b*x+self.c-distance_from_throat))
            x0 = array(float(self.throat_radius + (self.exit_radius-self.throat_radius)*(distance_from_throat/self.length)))
            radius = scipy.optimize.fsolve(func,x0)
        return radius

    @property
    def exit_area(self):
        return self.a_t * self.eps

    @property
    def throat_radius(self):
        return sqrt(self.a_t / pi)

    @property
    def exit_radius(self):
        return sqrt(self.exit_area / pi)


class CombustionChamber(Structure):
    def __init__(self, throat_area: float, combustion_chamber_pressure, convergent_half_angle: float,
                 in_radians: bool = True, throat_bend_ratio: float = 1., chamber_bend_ratio: float = .8,
                 propellant_mix: Optional[str] = None, area_ratio_cc_throat: Optional[float] = None,
                 characteristic_length: Optional[float] = None, **kwargs):
        # Initialize other parameters
        self.a_t = throat_area  # m2
        self.r_t = sqrt(self.a_t / pi)  # m
        self.p_cc = combustion_chamber_pressure  # Pa
        self.cc_bend = chamber_bend_ratio if chamber_bend_ratio is not None else 1
        self.th_bend = throat_bend_ratio if throat_bend_ratio is not None else .8
        if convergent_half_angle is not None:
            self.ha_conv = convergent_half_angle if in_radians else radians(convergent_half_angle)
        else:
            self.ha_conv = radians(30)
        if self.ha_conv > pi / 2 or self.ha_conv < 0:
            raise ValueError(
                f'The divergence half angle of the convergent must be between 0 and {pi / 2:.3f} ({self.ha_conv:.3f} given)')

        # Set defaults for optional parameters
        options = {
            'LOX/GH2': 0.635,
            'LOX/LH2': 0.89,
            'LOX/RP1': 1.145
        }
        if characteristic_length is None:
            try:
                self.l_star = options[propellant_mix]
            except KeyError:
                raise KeyError(
                    f'The specified propellant mix [{propellant_mix}] does not have a default characteristic length for the combustion chamber, specify one manually or select a propellant mix that does [{options.keys}]')
        else:
            self.l_star = characteristic_length

        if combustion_chamber_pressure is not None:
            self.ac_at = area_ratio_cc_throat
        else:
            # Humble 1995 p.222
            self.ac_at = (8.0 * (2 * self.r_t) ** 2.4 + 1.25)

        super().__init__(**kwargs)

    @property
    def volume(self):
        return self.l_star * self.a_t

    @property
    def mass(self):
        return self.md * self.sf / self.sy * 2 * self.volume * self.p_cc

    @property
    def length(self):
        return self.volume / self.area

    def convergent_radius(self, distance_from_throat):
        length_q = self.r_u * sin(self.ha_conv)
        radius_q = self.r_t + self.r_u * (1 - cos(self.ha_conv))
        radius_p = self.radius - self.r_a * (1 - cos(self.ha_conv))
        length_p = self.length_p
        if distance_from_throat > self.convergent_length or distance_from_throat < 0:
            raise ValueError(
                f'Radius of the convergent cannot be calculated above its length ({self.convergent_length:.4e} or downstream of the throat')
        elif distance_from_throat > length_p:
            distance = self.convergent_length - distance_from_throat
            alpha = asin(distance / self.r_a) / 2
            radius = self.radius - distance * tan(alpha)
        elif distance_from_throat > length_q:
            radius = radius_q + (distance_from_throat - length_q) * tan(self.ha_conv)
        else:
            alpha = asin(distance_from_throat / self.r_u) / 2
            radius = self.r_t + distance_from_throat * tan(alpha)
        return radius

    @property
    def area(self):
        return self.ac_at * self.a_t

    @property
    def radius(self):
        return sqrt(self.area / pi)

    @property
    def convergent_length(self):
        return self.length_p + self.r_a * sin(self.ha_conv)

    @property
    def length_p(self):
        z = self.radius - self.r_t - (self.r_u + self.r_a) * (1 - cos(self.ha_conv))
        if z < 0:
            raise ValueError(
                f'The length of the straight part of the convergent is calculated to be less than zero, please change the bend ratios or divergence half angle')
        return self.r_u * sin(self.ha_conv) + z / tan(self.ha_conv)

    @property
    def total_length(self):
        return self.convergent_length + self.length

    @property
    def r_u(self):
        # Longitudinal radius at throat [m]
        return self.th_bend * self.r_t

    @property
    def r_a(self):
        # Longitudinal radius between cylindrical chamber and convergent
        return self.cc_bend * self.radius


class ComplexCombustionChamber(CombustionChamber):
    def __init__(self, average_wall_temperature: Optional[float], **kwargs):
        self.t_w_avg = 700 if average_wall_temperature is None else average_wall_temperature  # Kelvin
        super().__init__(**kwargs)
        # Inconel 600 Data
        temperature = [
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
            1373
        ]
        ultimate_tensile_strength = [
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
            27
        ]
        sigma_ult_function = interp1d(temperature, ultimate_tensile_strength)
        self.sy = sigma_ult_function(self.t_w_avg)