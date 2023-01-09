from dataclasses import dataclass, field
from math import pi
from typing import Optional, Callable

import numpy
from scipy import constants as constants

from EngineComponents.Base.ThrustChamber import ThrustChamber
from EngineComponents.Abstract.FlowState import FlowState
from EngineFunctions.EmpiricalRelations import get_netto_average_wall_radiative_heat_flux, \
    get_hot_gas_convective_heat_transfer_coefficient


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
        return get_netto_average_wall_radiative_heat_flux(
            combustion_temperature=self.combustion_temperature,
            maximum_wall_temperature=self.maximum_wall_temperature,
            thrust_chamber_wall_emissivity=self.thrust_chamber_wall_emissivity,
            hot_gas_emissivity=self.hot_gas_emissivity
        )

    @property
    def theoretical_total_radiative_heat_transfer(self):  # [W]
        return self.netto_average_wall_radiative_heat_flux * self.thrust_chamber.surface_area

    @property
    def radiative_factor(self):
        return self.theoretical_total_radiative_heat_transfer / self.theoretical_total_convective_heat_transfer


@dataclass
class ConvectiveHeatTransfer:
    thrust_chamber: ThrustChamber
    combustion_chamber_flow_state: FlowState
    hot_gas_emissivity: float  # [-]
    maximum_wall_temperature: float  # [K]
    recovery_factor: float  # [-]
    post_injection_build_up_ratio: float = 0.25  # [-]
    # post_injection_build_up_ratio: Heat transfer to combustion chamber wall is assumed zero at injector face and
    # builds up to a constant heat transfer based on combustion temperature as reference temperature. At which point
    # this constant heat transfer is reached as percentage of total chamber length, is determined by
    # post_injection_build_up_ratio. Default value based on rough estimate from Perakis2021
    verbose: bool = True

    heat_transfer_func: Callable = field(init=False, repr=False)
    total_convective_heat_transfer: float = field(init=False, repr=False)
    _interpolation_num: float = 120

    def __post_init__(self):
        self.init_heat_transfer()

    @property
    def combustion_temperature(self):
        return self.combustion_chamber_flow_state.temperature

    @property
    def heat_capacity_ratio(self):
        return self.combustion_chamber_flow_state.heat_capacity_ratio

    @property
    def min_distance_from_throat(self):
        return self.thrust_chamber.min_distance_from_throat

    @property
    def max_distance_from_throat(self):
        return self.thrust_chamber.max_distance_from_throat

    def get_convective_heat_transfer_coefficient(self, distance_from_throat: float):
        return get_hot_gas_convective_heat_transfer_coefficient(
            mode="ModifiedBartz",
            mass_flow=self.combustion_chamber_flow_state.mass_flow,
            local_diameter=2 * self.thrust_chamber.get_radius(distance_from_throat),
            dynamic_viscosity=self.combustion_chamber_flow_state.dynamic_viscosity,
            specific_heat_capacity=self.combustion_chamber_flow_state.specific_heat_capacity,
            prandtl_number=self.combustion_chamber_flow_state.prandtl_number,
            film_temp=self.get_film_temperature(distance_from_throat),
            stagnation_temp=self.combustion_temperature,
        )

    def get_static_temp(self, distance_from_throat: float):
        m = self.thrust_chamber.get_mach(distance_from_throat)
        return self.combustion_temperature / (1 + (self.heat_capacity_ratio - 1) / 2 * m ** 2)

    def get_adiabatic_wall_temp(self, distance_from_throat: float):
        m = self.thrust_chamber.get_mach(distance_from_throat)
        y = self.heat_capacity_ratio
        r = self.recovery_factor
        factor1 = 1 + (y - 1) / 2 * m ** 2 * r
        factor2 = 1 + (y - 1) / 2 * m ** 2
        return self.combustion_temperature * factor1 / factor2

    def get_film_temperature(self, distance_from_throat: float):
        t_aw = self.get_adiabatic_wall_temp(distance_from_throat)
        t = self.get_static_temp(distance_from_throat)
        t_w = self.maximum_wall_temperature
        return float(.5 * t_w + .28 * t + t_aw * .22)

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
        self.heat_transfer_func = lambda x: numpy.interp(x, xs, ys)

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
