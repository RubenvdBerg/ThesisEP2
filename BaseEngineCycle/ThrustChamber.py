import warnings
from math import pi
from typing import Callable, Optional
import scipy.integrate
import scipy.optimize
from matplotlib import pyplot as plt
from numpy import array, isclose, linspace
from dataclasses import dataclass
from BaseEngineCycle.CombustionChamber import CombustionChamber, Injector
from BaseEngineCycle.Nozzle import Nozzle
from irt import get_expansion_ratio


@dataclass
class ThrustChamber:
    nozzle: Nozzle
    chamber: CombustionChamber
    injector: Injector
    heat_capacity_ratio: float  # [-]

    @property
    def min_distance_from_throat(self):
        return float((self.nozzle.conv_length + self.chamber.length) * -1)

    @property
    def max_distance_from_throat(self):
        return self.nozzle.div_length

    @property
    def length(self):
        return self.max_distance_from_throat - self.min_distance_from_throat

    @property
    def throat_distance_tuple(self):
        return self.min_distance_from_throat, self.max_distance_from_throat

    @property
    def throat_area(self):
        return self.nozzle.throat_area

    def get_radius(self, distance_from_throat):
        if distance_from_throat < self.min_distance_from_throat:
            raise ValueError(
                f'Radius of thrust chamber cannot be calculated before the injector plate ({self.min_distance_from_throat:.4e} m before the throat)')
        elif distance_from_throat < -self.nozzle.conv_length:
            return self.chamber.radius
        else:
            return self.nozzle.get_radius(distance_from_throat)

    @property
    def surface(self):
        result = scipy.integrate.quad(lambda x: self.get_radius(x) * 2 * pi, *self.throat_distance_tuple)
        return result[0]

    def get_distance_for_divergent_expansion_ratio(self, expansion_ratio):
        # Ugly way of getting distance from throat at a certain expansion ratio for the same nozzle
        nozzle_copy = deepcopy(self.nozzle)
        nozzle_copy.expansion_ratio = expansion_ratio
        return nozzle_copy.div_length

    def get_mach(self, distance_from_throat):
        radius = self.get_radius(distance_from_throat)
        expansion_ratio = pi * radius ** 2 / self.throat_area
        x0 = .1 if distance_from_throat < 0 else 10

        def func(mach_number):
            return array([float(get_expansion_ratio(mach_number, self.heat_capacity_ratio) - expansion_ratio)])

        mach = scipy.optimize.fsolve(func, array([float(x0)]))[0]
        check = isclose(func(mach), [0.0])
        if not check:
            warnings.warn('Implicitly calculated Mach number is not within tolerances')
        return float(mach)

    def show_mach(self, **kwargs):
        self.distance_plot(self.get_mach, 'Mach [-]', **kwargs)

    def show_contour(self, **kwargs):
        self.distance_plot(self.get_radius, 'Radius [m]', **kwargs)

    def distance_plot(self, func: Callable, ylabel: str, num=300, ytick_function: Optional[Callable] = None, distance_tuple: Optional[tuple] = None):
        if distance_tuple is None:
            distance_tuple = self.throat_distance_tuple
        distances = list(linspace(*distance_tuple, num))
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
