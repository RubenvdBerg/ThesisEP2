from dataclasses import dataclass
from functools import cache
from math import log, pi
from typing import Optional
import scipy.integrate
from scipy import constants as constants
from BaseEngineCycle.ThrustChamber import ThrustChamber


@dataclass
class Coolant:
    heat_capacity_liquid: float  # [J/mol]
    heat_capacity_gas: float  # [J/mol]
    heat_of_vaporization: float  # [J/mol]
    molar_mass: float  # [kg/mol]!!
    boiling_temperature_1_bar: float  # [K]

    # Coolant properties
    @property
    @cache
    def specific_heat_capacity_liquid(self):
        return self.heat_capacity_liquid / self.molar_mass

    @property
    @cache
    def specific_heat_capacity_gas(self):
        return self.heat_capacity_gas / self.molar_mass

    @property
    @cache
    def specific_heat_of_vaporization(self):
        return self.heat_of_vaporization / self.molar_mass

    @cache
    def boiling_temperature(self, pressure):
        # Clausius-Clapeyron estimation
        T0 = self.boiling_temperature_1_bar
        Dh = self.cooolant_heat_of_vaporization
        return (T0 ** -1 - constants.R / Dh * log(pressure)) ** -1

    @cache
    def start_boiling_enthalpy(self, pressure):
        return self.boiling_temperature(pressure) * self.heat_capacity_liquid

    @cache
    def end_boiling_enthalpy(self, pressure):
        return self.start_boiling_enthalpy(pressure) + self.heat_of_vaporization


@dataclass
class CoolingChannels:
    coolant: Coolant
    total_heat_transfer: float  # [W]
    outlet_pressure: float  # [Pa]
    inlet_temperature: float  # [K]
    mass_flow: float  # [kg/s]
    pressure_drop: Optional[float] = None  # [Pa]
    _pressure_drop_ratio: float = .15  # [-]

    @property
    @cache
    def inlet_pressure(self):
        if self.pressure_drop is None:
            # Humble 1995 p.209 suggest pressure drop to be 10% - 20% of chamber/outlet pressure
            return self.outlet_pressure * (1 + self._pressure_drop_ratio)
        else:
            return self.outlet_pressure - self.pressure_drop

    @property
    @cache
    def inlet_enthalpy(self):
        if self.coolant.boiling_temperature(self.inlet_pressure) < self.inlet_temperature:
            raise ValueError(
                'The boiling temperature of the coolant is below the inlet temperature, liquid cooling expected')
        return self.inlet_temperature * self.coolant.specific_heat_capacity_liquid

    @property
    @cache
    def enthalpy_increase(self):
        return self.total_heat_transfer / self.mass_flow

    @property
    @cache
    def outlet_enthalpy(self):
        return self.inlet_enthalpy + self.enthalpy_increase

    @property
    @cache
    def outlet_temperature(self):
        cp_g = self.coolant.specific_heat_capacity_gas
        cp_l = self.coolant.specific_heat_capacity_liquid
        temp_boil = self.coolant.boiling_temperature(self.inlet_pressure)
        h_vap_start = self.coolant.start_boiling_enthalpy(self.inlet_pressure)
        h_vap_end = self.coolant.end_boiling_enthalpy(self.inlet_pressure)

        if self.outlet_enthalpy > h_vap_end:
            return (self.outlet_enthalpy - h_vap_end) / cp_g + temp_boil
        elif self.outlet_enthalpy > h_vap_start:
            return temp_boil
        else:
            return self.outlet_enthalpy / cp_l

    @property
    @cache
    def coolant_outlet_temperature(self):
        return self.inlet_temperature + self.temperature_increase


@dataclass
class HeatExchanger:
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
    prandtl_number: Optional[float] = None  # [-]
    recovery_factor: Optional[float] = None  # [-]

    def __post_init__(self):
        # Setting recovery factor with turbulent estimate if not provided
        self.recovery_factor = self.turbulent_recovery_factor if self.recovery_factor is None else self.recovery_factor
        self.prandtl_number = self.prandtl_number_estimate if self.prandtl_number is None else self.prandtl_number
    
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
    def total_radiative_heat_transfer(self):
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
            if distance_from_throat < -self.thrust_chamber.chamber.conv_length:
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
        if distance_from_throat < -self.thrust_chamber.nozzle.conv_length:
            temp_ref = self.combustion_temperature
        else:
            temp_ref = self.get_adiabatic_wall_temp(distance_from_throat)
        return coefficient * (temp_ref - self.maximum_wall_temperature)

    @property
    def total_convective_heat_transfer(self):
        result = scipy.integrate.quad(
            lambda x: self.get_convective_heat_flux(x) * self.thrust_chamber.get_radius(x) * 2 * pi,
            *self.thrust_chamber.throat_distance_tuple)
        # if self.verbose:
        #     print(
        #         f'Total Convective Heat Transfer estimated with a estimated error of {result[1] / result[0] * 100:.8f}%')
        return float(result[0])

    @property
    def total_heat_transfer(self):
        return self.total_convective_heat_transfer + self.total_radiative_heat_transfer

    def show_heat_flux_coefficient(self, **kwargs):
        self.thrust_chamber.distance_plot(func=self.get_convective_heat_transfer_coefficient,
                                          ylabel=r'Convective Heat Transfer Coefficient [$kW$/$(m^2K)$]',
                                          ytick_function=lambda x: f'{x * 1e-3:.0f}',
                                          **kwargs)

    def show_heat_flux(self, **kwargs):
        self.thrust_chamber.distance_plot(func=self.get_convective_heat_flux,
                                          ylabel=r'Convective Heat Flux [$MW$/$m^2$]',
                                          ytick_function=lambda x: f'{x * 1e-6:.0f}',
                                          **kwargs)

    def show_adiabatic_wall_temp(self, **kwargs):
        self.thrust_chamber.distance_plot(func=self.get_adiabatic_wall_temp,
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