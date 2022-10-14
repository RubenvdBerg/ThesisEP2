import warnings
from math import pi, sqrt
from typing import Callable, Optional
import scipy.integrate
import scipy.optimize
from matplotlib import pyplot as plt
from numpy import array, isclose, linspace
from dataclasses import dataclass
from BaseEngineCycle.CombustionChamber import CombustionChamber
from BaseEngineCycle.Injector import Injector
from BaseEngineCycle.Nozzle import Nozzle
from irt import get_expansion_ratio, get_local_mach
from copy import deepcopy


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
        larger = distance_from_throat > self.max_distance_from_throat
        smaller = distance_from_throat < self.min_distance_from_throat
        if smaller or larger:
            raise ValueError(f'Distance from throat larger or smaller than {type(self)} bounds')
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
        if distance_from_throat < -self.nozzle.conv_length:
            return 0
        else:
            radius = self.get_radius(distance_from_throat)
            local_area_ratio = pi * radius ** 2 / self.throat_area
            is_subsonic = True if distance_from_throat < 0 else False
            return get_local_mach(local_area_ratio=local_area_ratio,
                                  is_subsonic=is_subsonic,
                                  heat_capacity_ratio=self.heat_capacity_ratio)

    def show_mach(self, **kwargs):
        self.distance_plot(self.get_mach, 'Mach [-]', **kwargs)

    def show_contour(self, **kwargs):
        self.distance_plot(self.get_radius, 'Radius [m]', **kwargs)

    def distance_plot(self, func: Callable, ylabel: str, num=300, ytick_function: Optional[Callable] = None,
                      distance_tuple: Optional[tuple] = None, distance_from_injector: bool = False):
        if distance_tuple is None:
            distance_tuple = self.throat_distance_tuple
        distances = list(linspace(*distance_tuple, num))
        values = [func(distance) for distance in distances]
        fig, ax = plt.subplots()
        if distance_from_injector:
            distances = [x - distances[0] for x in distances]
            xlabel = 'Distance from injector face [m]'
        else:
            xlabel = 'Distance from throat [m]'

        ax.plot(distances, values)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        if ytick_function is not None:
            ticks = ax.get_yticks().tolist()
            ax.set_yticks(ticks)
            ax.set_yticklabels([ytick_function(x) for x in ticks])
        plt.show()


@dataclass
class ThrustChamberSection(ThrustChamber):
    """Provide a section of a thrust chamber.

    Section bounds are determined by given distances from throat or by expansion ratios, which are translated to
    distance from the throat in the DIVERGENT! section of the nozzle.
    """
    min_distance: Optional[float] = None  # [m]
    max_distance: Optional[float] = None  # [m]
    min_distance_expansion_ratio: Optional[float] = None  # [-]
    max_distance_expansion_ratio: Optional[float] = None  # [-]

    def __post_init__(self):
        for minmax in ['min', 'max']:
            eps = getattr(self, f'{minmax}_distance_expansion_ratio')
            if eps is not None:
                distance = self.get_distance_for_divergent_expansion_ratio(eps)
                setattr(self, f'{minmax}_distance', distance)
                warnings.warn(
                    f'{minmax}_distance_expansion_ratio was provided, {minmax}_distance input will be ignored')

    @property
    def min_distance_from_throat(self):
        if self.min_distance is None:
            return super().min_distance_from_throat
        else:
            return self.min_distance

    @property
    def max_distance_from_throat(self):
        if self.max_distance is None:
            return super().max_distance_from_throat
        else:
            return self.max_distance
