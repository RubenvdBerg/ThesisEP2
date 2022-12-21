from dataclasses import dataclass, field
from math import pi
from typing import Optional, Callable

import numpy
import scipy.integrate
import scipy.interpolate
from scipy import constants as constants

from EngineCycles.BaseEngineCycle.ThrustChamber import ThrustChamber


@dataclass
class RadiativeHeatTransfer:
    thrust_chamber: ThrustChamber
    combustion_temperature: float  # [K}
    hot_gas_emissivity: float  # [-]
    maximum_wall_temperature: float  # [K]
    thrust_chamber_wall_emissivity: float  # [-]
    theoretical_total_convective_heat_transfer: float  # [-] Theoretical total heat transfer to complete surface of the thrust chamber

    @property
    def netto_average_wall_radiative_heat_flux(self):  # q_rad [W/m2]
        # Heat Transfer Handbook, A. Bejan 2003, Eq. 8.69
        tc = self.combustion_temperature
        tw = self.maximum_wall_temperature
        e_cw = self.thrust_chamber_wall_emissivity
        e_hg = self.hot_gas_emissivity
        return constants.sigma * (tc ** 4 - tw ** 4) / (1 / e_hg + (1 / e_cw) - 1)

    @property
    def theoretical_total_radiative_heat_transfer(self):  # [W]
        return self.netto_average_wall_radiative_heat_flux * self.thrust_chamber.surface_area

    @property
    def radiative_factor(self):
        return self.theoretical_total_radiative_heat_transfer / self.theoretical_total_convective_heat_transfer


@dataclass
class ConvectiveHeatTransfer:
    thrust_chamber: ThrustChamber
    # Properties of hot gas in combustion chamber
    combustion_temperature: float  # [K}
    combustion_chamber_pressure: float  # [Pa]
    mass_flow: float  # [kg/s]
    dynamic_viscosity: float  # [Pa*s]
    specific_heat_capacity: float  # [J/(kg*K)]
    hot_gas_emissivity: float  # [-]
    heat_capacity_ratio: float  # [-]
    maximum_wall_temperature: float  # [K]
    convective_coefficient_mode: str
    post_injection_build_up_ratio: Optional[float] = None  # [-]
    prandtl_number: Optional[float] = None  # [-]
    recovery_factor: Optional[float] = None  # [-]
    verbose: bool = True

    heat_transfer_func: Callable = field(init=False, repr=False)
    total_convective_heat_transfer: float = field(init=False, repr=False)
    _interpolation_num: float = 120

    def __post_init__(self):
        # Setting recovery factor with turbulent estimate if not provided
        # ORDER OF THESE LINES IS IMPORTANT
        if self.prandtl_number is None:
            self.prandtl_number = self.prandtl_number_estimate
        if self.recovery_factor is None:
            self.recovery_factor = self.turbulent_recovery_factor
        if self.post_injection_build_up_ratio is None:
            # Heat transfer to combustion chamber wall is assumed zero at injector face and builds up to a constant
            # heat transfer based on combustion temperature as reference temperature. At which point this constant heat
            # transfer is reached as percentage of total chamber length, is determined by post_injection_build_up_ratio
            # Rough estimate based on Perakis2021
            self.post_injection_build_up_ratio = 0.25

        self.init_heat_transfer()

    @property
    def min_distance_from_throat(self):
        return self.thrust_chamber.min_distance_from_throat

    @property
    def max_distance_from_throat(self):
        return self.thrust_chamber.max_distance_from_throat

    @property
    def turbulent_recovery_factor(self):
        # Zandbergen 2017 p.160
        return self.prandtl_number ** (1 / 3)

    @property
    def prandtl_number_estimate(self):
        # Zandbergen 2017 p.159
        y = self.heat_capacity_ratio
        return 4 * y / (9 * y - 5)

    @staticmethod
    def convective_heat_transfer_coefficient(mode: str, mass_flow: float, diameter: float, dynamic_viscosity: float,
                                             specific_heat_capacity: float, prandtl_number: float,
                                             total_temperature: Optional[float] = None,
                                             film_temperature: Optional[float] = None) -> float:
        mf = mass_flow
        di = diameter
        mu = dynamic_viscosity
        cp = specific_heat_capacity
        pr = prandtl_number
        t0 = total_temperature
        tf = film_temperature
        # Zandbergen 2017 p.161
        modes = ["Modified Bartz", "Cornelisse", "CornelisseNozzle", "Standard Bartz"]
        if mode == modes[0]:
            return 0.026 * 1.213 * mf ** .8 * di ** -1.8 * mu ** .2 * cp * pr ** -.6 * (t0 / tf) ** .68
        elif mode == modes[1]:
            return 0.023 * 1.213 * mf ** .8 * di ** -1.8 * mu ** .2 * cp * pr ** float(-2 / 3)
        elif mode == modes[2]:
            return 0.026 * 1.213 * mf ** .8 * di ** -1.8 * mu ** .2 * cp * pr ** float(-2 / 3) * (t0 / tf) ** .68
        elif mode == modes[3]:
            raise NotImplementedError(
                "Convective heat transfer for the standard bartz equation has not been implemented")
        else:
            raise ValueError(
                f"Improper convective_mode given for calculation of the convective heat transfer, pick one of {modes}")

    def get_convective_heat_transfer_coefficient(self, distance_from_throat: float):
        diameter = 2 * self.thrust_chamber.get_radius(distance_from_throat)

        # Cornelisse is not suited for the nozzle without a temperature correction
        if self.convective_coefficient_mode == "Cornelisse":
            if distance_from_throat < -self.thrust_chamber.nozzle.conv_length:
                mode = "Cornelisse"
            else:
                mode = "CornelisseNozzle"
        elif self.convective_coefficient_mode == "Modified Bartz":
            mode = self.convective_coefficient_mode
        else:
            mode = None
        return self.convective_heat_transfer_coefficient(
            mode=mode,
            mass_flow=self.mass_flow,
            diameter=diameter,
            dynamic_viscosity=self.dynamic_viscosity,
            specific_heat_capacity=self.specific_heat_capacity,
            prandtl_number=self.prandtl_number,
            film_temperature=self.get_film_temperature(distance_from_throat),
            total_temperature=self.combustion_temperature
        )

    def get_static_temp(self, distance_from_throat: float):
        m = self.thrust_chamber.get_mach(distance_from_throat)
        return self.combustion_temperature / (1 + (self.heat_capacity_ratio - 1) / 2 * m ** 2)

    def get_film_temperature(self, distance_from_throat: float):
        t_aw = self.get_adiabatic_wall_temp(distance_from_throat)
        t = self.get_static_temp(distance_from_throat)
        t_w = self.maximum_wall_temperature
        return float(.5 * t_w + .28 * t + t_aw * .22)

    def get_adiabatic_wall_temp(self, distance_from_throat: float):
        m = self.thrust_chamber.get_mach(distance_from_throat)
        y = self.heat_capacity_ratio
        r = self.recovery_factor
        factor1 = 1 + (y - 1) / 2 * m ** 2 * r
        factor2 = 1 + (y - 1) / 2 * m ** 2
        return self.combustion_temperature * factor1 / factor2

    def get_convective_heat_flux(self, distance_from_throat: float):
        coefficient = self.get_convective_heat_transfer_coefficient(distance_from_throat)
        len_cc = self.thrust_chamber.chamber.length
        dist_min = self.thrust_chamber.min_distance_from_throat
        len_conv = self.thrust_chamber.nozzle.conv_length
        r_build_up = self.post_injection_build_up_ratio
        if distance_from_throat < -len_conv:
            temp_ref = self.combustion_temperature
            # Distance from injector divided by total chamber length
            inj_distance_ratio = (distance_from_throat - dist_min) / len_cc
            if inj_distance_ratio < r_build_up:
                fact_distance = ((distance_from_throat - dist_min) / (r_build_up * len_cc))
            else:
                fact_distance = 1
            temp_eff = (temp_ref - self.maximum_wall_temperature) * fact_distance
        else:
            temp_ref = self.get_adiabatic_wall_temp(distance_from_throat)
            temp_eff = temp_ref - self.maximum_wall_temperature
        return coefficient * temp_eff

    def get_convective_heat_transfer_per_axial_meter(self, distance_from_throat: float):
        return (2 * pi * self.get_convective_heat_flux(distance_from_throat)
                * self.thrust_chamber.get_radius(distance_from_throat))

    @property
    def distance_tuple(self):
        return self.min_distance_from_throat, self.max_distance_from_throat

    def init_heat_transfer(self):
        xs = numpy.linspace(*self.distance_tuple, self._interpolation_num)
        dx = (self.max_distance_from_throat - self.min_distance_from_throat) / (self._interpolation_num - 1)
        ys = numpy.array([self.get_convective_heat_transfer_per_axial_meter(x) for x in xs])
        self.total_convective_heat_transfer = sum(dx * ys)
        self.heat_transfer_func = lambda x : numpy.interp(x, xs, ys)

    def distance_plot(self, **kwargs):
        self.thrust_chamber.distance_plot(**kwargs, distance_tuple=self.distance_tuple)

    def show_heat_flux_coefficient(self, **kwargs):
        self.distance_plot(func=self.get_convective_heat_transfer_coefficient,
                           ylabel=r'Convective Heat Transfer Coefficient [$kW$/$(m^2K)$]',
                           ytick_function=lambda x: f'{x * 1e-3:.0f}',
                           **kwargs)

    def show_heat_flux(self, **kwargs):
        self.distance_plot(func=self.get_convective_heat_flux,
                           ylabel=r'Convective Heat Flux [$MW$/$m^2$]',
                           ytick_function=lambda x: f'{x * 1e-6:.0f}',
                           **kwargs)

    def show_heat_transfer(self, **kwargs):
        self.distance_plot(func=self.get_convective_heat_transfer_per_axial_meter,
                           ylabel=r'Convective Heat Transfer/m [$MW$/$m$]',
                           ytick_function=lambda x: f'{x * 1e-6:.1f}',
                           **kwargs)

    def show_adiabatic_wall_temp(self, **kwargs):
        self.distance_plot(func=self.get_adiabatic_wall_temp,
                           ylabel=r'Adiabatic Wall Temperature [$K$]',
                           **kwargs)


@dataclass
class HeatTransferSection(ConvectiveHeatTransfer):
    min_distance_section: Optional[float] = None  # [m]
    max_distance_section: Optional[float] = None  # [m]
    radiative_heat_transfer_factor: float = 0

    @property
    def min_distance_from_throat(self):
        if self.min_distance_section is None:
            return self.thrust_chamber.min_distance_from_throat
        else:
            return self.min_distance_section

    @property
    def max_distance_from_throat(self):
        if self.max_distance_section is None:
            return self.thrust_chamber.max_distance_from_throat
        else:
            return self.max_distance_section

    def total_heat_flux(self, distance_from_throat: float):
        return self.get_convective_heat_flux(distance_from_throat) * (1 + self.radiative_heat_transfer_factor)

    @property
    def total_heat_transfer(self):
        return self.total_convective_heat_transfer * (self.radiative_heat_transfer_factor + 1)