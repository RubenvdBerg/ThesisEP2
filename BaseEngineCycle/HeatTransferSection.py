from dataclasses import dataclass
from math import pi
from typing import Optional
from copy import deepcopy

import scipy.integrate
from scipy import constants as constants

from BaseEngineCycle.ThrustChamber import ThrustChamber
from BaseEngineCycle.Nozzle import Nozzle


@dataclass
class HeatTransferSection:
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

    thrust_chamber_wall_emissivity: float  # [-]
    convective_coefficient_mode: str
    min_distance_section: Optional[float] = None  # [m]
    max_distance_section: Optional[float] = None  # [m]
    post_injection_build_up_ratio: Optional[float] = None  # [-]
    prandtl_number: Optional[float] = None  # [-]
    recovery_factor: Optional[float] = None  # [-]
    verbose: bool = True

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

    @property
    def laminar_recovery_factor(self):
        return self.prandtl_number ** .5

    @property
    def turbulent_recovery_factor(self):
        # Zandbergen 2017 p.160
        return self.prandtl_number ** (1 / 3)

    @property
    def prandtl_number_estimate(self):
        # Zandbergen 2017 p.159
        y = self.heat_capacity_ratio
        return 4 * y / (9 * y - 5)

    @property
    def netto_average_wall_radiative_heat_flux(self):  # q_rad [W/m2]
        # Heat Transfer Handbook, A. Bejan 2003, Eq. 8.69
        tc = self.combustion_temperature
        tw = self.maximum_wall_temperature
        e_cw = self.thrust_chamber_wall_emissivity
        e_hg = self.hot_gas_emissivity
        return constants.sigma * (tc ** 4 - tw ** 4) / (1 / e_hg + (1 / e_cw) - 1)

    @property
    def total_radiative_heat_transfer(self):  # [W]
        return self.netto_average_wall_radiative_heat_flux * self.thrust_chamber.surface

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
        t_0 = self.get_static_temp(distance_from_throat)
        return float((t_0 + self.maximum_wall_temperature) / 2)

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

    @property
    def distance_tuple(self):
        return self.min_distance_from_throat, self.max_distance_from_throat

    @property
    def total_convective_heat_transfer(self):  # [W]
        result = scipy.integrate.quad(
            lambda x: self.get_convective_heat_flux(x) * self.thrust_chamber.get_radius(x) * 2 * pi,
            *self.distance_tuple)
        # if self.verbose:
        #     print(
        #         f'Total Convective Heat Transfer estimated with a estimated error of {result[1] / result[0] * 100:.8f}%')
        return float(result[0])

    @property
    def total_heat_transfer(self):  # [W]
        return self.total_convective_heat_transfer + self.total_radiative_heat_transfer

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


    def show_adiabatic_wall_temp(self, **kwargs):
        self.distance_plot(func=self.get_adiabatic_wall_temp,
                           ylabel=r'Adiabatic Wall Temperature [$K$]',
                           **kwargs)

    # def convective_heat_flux(self):
    #     if self.convective_mode == "Modified Bartz":
    #         return 0.026 * 1.213 * self.mass_flow**.8 * self.diameter**-1.8 * self.mu**.2 * self.cp * self.prandtl**-.6 * (self.tc / self.tf)**.68
    #     elif self.convective_mode == "Cornellisse":
    #         return 0.023 * 1.213 * self.mass_flow**.8 * self.diameter**-1.8 * self.mu**.2 * self.cp * self.prandtl**-2/3
    #     elif self.convective_mode == "Standard Bartz":
    #         raise NotImplementedError("Convective heat transfer for the standard bartz equation has not been implemented")
    #     else:
    #         raise ValueError("Improper convective_mode given for calculation of the convective heat transfer")
