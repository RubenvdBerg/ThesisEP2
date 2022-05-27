from math import pi, log, radians, sqrt
from functools import cache
import scipy.optimize
import scipy.constants as constants
import scipy.integrate as integrate
from BaseEngineCycle.ThrustChamber import ThrustChamber, Nozzle, CombustionChamber
from cea import get_cea_values
from typing import Optional, Literal
from irt import get_kerckhove, get_expansion_ratio
from dataclasses import dataclass
import warnings


@dataclass
class Structure:
    material_density: float  # [kg/m3]
    safety_factor: float  # [-]
    yield_strength: float  # [Pa]


@dataclass
class Tank(Structure):
    max_acceleration: float  # m/s2
    ullage_factor: float  # -
    initial_pressure: float  # Pa
    propellant: Propellant
    kwak_fix_cycle_type: str
    pressurant_tank: Optional[PressurantTank] = None
    kwak_fix: bool = False

    @property
    def radius(self):
        return (3 * self.volume / (4 * pi)) ** (1 / 3)

    @property
    def initial_head(self):
        unused_volume = (self.ullage_factor - 1) * self.propellant.volume
        rest_height = (unused_volume * 3 / (2 * pi)) ** (1 / 3)
        if self.kwak_fix:
            if self.kwak_fix_cycle_type == 'ep':
                if self.propellant.type == 'oxidizer':
                    return 1.91839449096392
                elif self.propellant.type == 'fuel':
                    return 1.65560478870526
            elif self.kwak_fix_cycle_type == 'gg':
                if self.propellant.type == 'oxidizer':
                    return 1.92539045861846
                elif self.propellant.type == 'fuel':
                    return 1.71033897378923
        return 2 * self.radius - rest_height

    @property
    def total_lower_pressure(self):
        return self.initial_pressure + self.propellant.density * self.max_acceleration * self.initial_head

    @property
    def total_upper_pressure(self):
        return self.initial_pressure + self.propellant.density * self.max_acceleration * (self.initial_head -
                                                                                          self.radius)

    @property
    def mass(self):  # Spherical Tanks
        return (self.safety_factor * 3 * self.material_density / (4 * self.yield_strength) *
                self.volume * (self.total_upper_pressure + self.total_lower_pressure))

    @property
    def volume(self):
        if self.propellant.type == 'oxidizer':
            if self.pressurant_tank.volume is not None:
                return self.propellant.volume * self.ullage_factor + self.pressurant_tank.volume
            else:
                warnings.warn(
                    'No pressurant tank volume was given. Assumed pressurant tank is not submerged in oxygen tank')
                return self.propellant.volume * self.ullage_factor
        elif self.propellant.type == 'fuel':
            return self.propellant.volume * self.ullage_factor


@dataclass
class Propellant:
    type: Literal['oxidizer', 'fuel']
    mass_flow: float  # [kg/s]
    burn_time: float  # [s]
    density: float  # [kg/m3]
    margin_factor: float  # [-]

    @property
    def mass(self):
        return self.mass_flow * self.burn_time * self.margin_factor

    @property
    def volume(self):
        return self.mass / self.density

    @property
    def volumetric_flow(self):
        return self.mass_flow / self.density


@dataclass
class Pressurant:
    fuel: Propellant
    oxidizer: Propellant
    fuel_tank: Tank
    oxidizer_tank: Tank
    margin_factor: float  # [-]
    structural_factor: float  # [-]
    initial_pressure: float  # [Pa]
    final_pressure: float  # [Pa]
    heat_capacity_ratio: float  # [-]
    molar_mass: float  # [kg/mol]!!
    initial_temperature: float  # [K]

    @property
    def specific_gas_constant(self):
        return constants.R / self.molar_mass

    @property
    def mass(self):
        fact_m = self.margin_factor
        fact_u = self.fuel_tank.ullage_factor
        y = self.heat_capacity_ratio
        R_sp = self.specific_gas_constant
        T_0 = self.initial_temperature
        otp = self.oxidizer_tank.initial_pressure
        ftp = self.fuel_tank.initial_pressure
        ov = self.oxidizer.volume
        fv = self.fuel_volume
        p1 = self.initial_pressure
        p0 = self.final_pressure
        return fact_m * fact_u * y / (R_sp * T_0) * (otp * ov + ftp * fv) / (1 - (p1 / p0))


@dataclass
class PressurantTank(Structure):
    pressurant: Pressurant

    @property
    def initial_pressure(self):
        return self.pressurant.initial_pressure

    @property
    def final_pressure(self):
        return self.pressurant.final_pressure

    @property
    def volume(self):
        return (self.pressurant.mass * self.pressurant.specific_gas_constant
                * self.pressurant.initial_temperature / self.pressurant.initial_pressure)

    @property
    def mass(self):
        return (self.safety_factor * 3 / 2 * self.volume * self.initial_pressure
                * self.material_density / self.yield_strength)


@dataclass
class Pump:
    propellant: Propellant
    mass_flow: float  # [kg/s]
    pressure_increase: float  # [Pa]
    efficiency: float  # [-]
    specific_power: float  # [W/kg]

    @property
    def volumetric_flow_rate(self):
        return self.mass_flow / self.propellant.density

    @property
    def power_required(self):
        return self.volumetric_flow_rate * self.pressure_increase / self.efficiency

    @property
    def mass(self):
        return self.power_required / self.specific_power


@dataclass
class Injector(Structure):
    combustion_chamber: CombustionChamber
    propellant_is_gas: bool

    @property
    def pressure_drop(self):
        # Mota 2008 -> Kesaev and Almeida 2005
        f = .4 if self.propellant_is_gas else .8
        return f * 10E2 * sqrt(10 * self.combustion_chamber.p_cc)

    @property
    def thickness(self):
        # Zandbergen -> 4 times pressure drop, ASME code for welded flat caps on internal pressure vessels
        diameter = 2 * sqrt(self.combustion_chamber.area / pi)
        pressure = 4 * self.pressure_drop
        c = 6 * 3.3 / 64
        return diameter * sqrt(c * pressure / self.yield_strength) * self.safety_factor

    @property
    def mass(self):
        return self.combustion_chamber.area * self.thickness * self.material_density


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
        if self.pressure_drop is not None:
            return self.outlet_pressure - self.pressure_drop
        else:
            # Humble 1995 p.209 suggest pressure drop to be 10% - 20% of chamber/outlet pressure
            return self.outlet_pressure * (1 + self._pressure_drop_ratio)

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
    cooling_channels: CoolingChannels
    # Properties of hot gas in combustion chamber
    combustion_temperature: float  # [K}
    combustion_chamber_pressure: float  # [Pa]
    mass_flow: float  # [kg/s]
    dynamic_viscosity: float  # [Pa*s]
    specific_heat_capacity: float  # [J/(kg*K)]
    prandtl_number: float  # [-]
    hot_gas_emissivity: float  # [-]
    heat_capacity_ratio: float  # [-]

    maximum_wall_temperature: float  # [K]
    thrust_chamber_wall_emissivity: float  # [-]
    convective_coefficient_mode: str
    recovery_factor: Optional[float] = None  # [-]

    @property
    def netto_average_wall_radiative_heat_flux(self):  # q_rad [W/m2]
        # Heat Transfer Handbook, A. Bejan 2003, Eq. 8.69
        tc = self.combustion_temperature
        tw = self.maximum_wall_temperature
        e_cw = self.thrust_chamber_wall_emissivity
        e_hg = self.hot_gas_emissivity
        return constants.sigma(tc ** 4 - tw ** 4) / (1 / e_hg + (1 / e_cw) - 1)

    @property
    def total_radiative_heat_transfer(self):
        return self.netto_average_wall_radiative_heat_flux * self.thrust_chamber.surface

    @property
    def laminar_recovery_factor(self):
        return self.prandtl_number ** .5

    @property
    def turbulent_recovery_factor(self):
        # Zandbergen 2017 p.160
        return self.prandtl_number ** (1 / 3)

    def get_convective_heat_transfer_coefficient(self, distance_from_throat: float):
        diameter = 2 * self.thrust_chamber.get_radius(distance_from_throat)
        if self.convective_coefficient_mode == "Cornelisse":
            if distance_from_throat < -self.thrust_chamber.cc.convergent_length:
                mode = "Cornelisse"
            else:
                mode = "CornelisseNozzle"
        elif self.convective_coefficient_mode == "Modified Bartz":
            mode = self.convective_coefficient_mode
        else:
            mode = None
        return convective_heat_transfer_coefficient(
            mode=mode,
            mass_flow=self.mass_flow,
            diameter=diameter,
            dynamic_viscosity=self.dynamic_viscosity,
            specific_heat_capacity=self.specific_heat_capacity,
            prandtl_number=self.prandtl_number,
            film_temperature=self.get_film_temperature(distance_from_throat),
            total_temperature=self.combustion_temperature
        )

    def get_convective_heat_flux(self, distance_from_throat: float):
        coefficient = self.get_convective_heat_transfer_coefficient(distance_from_throat)
        if distance_from_throat < -self.thrust_chamber.cc.convergent_length:
            temp_ref = self.combustion_temperature
        else:
            temp_ref = self.get_adiabatic_wall_temp(distance_from_throat)
        return coefficient * (temp_ref - self.maximum_wall_temperature)

    def get_static_temp(self, distance_from_throat: float):
        m = self.thrust_chamber.get_mach(distance_from_throat)
        return self.tc / (1 + (self.heat_capacity_ratio - 1) / 2 * m ** 2)

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
    def total_heat_flux(self):
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

    @staticmethod
    def convective_heat_transfer_coefficient(mode: str, mass_flow: float, diameter: float, dynamic_viscosity: float,
                                             specific_heat_capacity: float, prandtl_number: float,
                                             total_temperature: float = None,
                                             film_temperature: float = None):
        # Zandbergen 2017 p.161
        modes = ["Modified Bartz", "Cornelisse", "CornelisseNozzle", "Standard Bartz"]
        if mode == modes[0]:
            return 0.026 * 1.213 * mass_flow ** .8 * diameter ** -1.8 * dynamic_viscosity ** .2 * specific_heat_capacity * prandtl_number ** -.6 * (
                    total_temperature / film_temperature) ** .68
        elif mode == modes[1]:
            return 0.023 * 1.213 * mass_flow ** .8 * diameter ** -1.8 * dynamic_viscosity ** .2 * specific_heat_capacity * prandtl_number ** float(
                -2 / 3)
        elif mode == modes[2]:
            return 0.026 * 1.213 * mass_flow ** .8 * diameter ** -1.8 * dynamic_viscosity ** .2 * specific_heat_capacity * prandtl_number ** float(
                -2 / 3) * (
                           total_temperature / film_temperature) ** .68
        elif mode == modes[3]:
            raise NotImplementedError(
                "Convective heat transfer for the standard bartz equation has not been implemented")
        else:
            raise ValueError(
                f"Improper convective_mode given for calculation of the convective heat transfer, pick one of {modes}")

    # def convective_heat_flux(self):
    #     if self.convective_mode == "Modified Bartz":
    #         return 0.026 * 1.213 * self.mass_flow**.8 * self.diameter**-1.8 * self.mu**.2 * self.cp * self.prandtl**-.6 * (self.tc / self.tf)**.68
    #     elif self.convective_mode == "Cornellisse":
    #         return 0.023 * 1.213 * self.mass_flow**.8 * self.diameter**-1.8 * self.mu**.2 * self.cp * self.prandtl**-2/3
    #     elif self.convective_mode == "Standard Bartz":
    #         raise NotImplementedError("Convective heat transfer for the standard bartz equation has not been implemented")
    #     else:
    #         raise ValueError("Improper convective_mode given for calculation of the convective heat transfer")


# def get_static_temp(combustion_temperature, heat_capacity_ratio, mach_number):
#     return combustion_temperature / (1 + (heat_capacity_ratio - 1) / 2 * mach_number ** 2)
#
#
# def get_film_temperature(static_temperature, wall_temperature):
#     return (static_temperature + wall_temperature) / 2
#
#
# def get_adiabatic_wall_temp(static_temp, heat_capacity_ratio, mach_number, recovery_factor):
#     return static_temp * (1 + (heat_capacity_ratio - 1) / 2 * mach_number ** 2) * recovery_factor
#
#
# def get_convective_heat_flux(heat_transfer_coefficient, reference_temperature, wall_temperature):
#     return heat_transfer_coefficient * (reference_temperature - wall_temperature)


@dataclass
class BaseEngineInput:
    thrust: float
    burn_time: float
    combustion_chamber_pressure: float
    oxidizer_name: str
    fuel_name: str
    is_frozen: bool
    exit_pressure: float
    max_acceleration: float
    heat_ratio_pressurant: float
    mass_mixture_ratio: float
    pressurant_initial_pressure: float
    pressurant_final_pressure: float
    oxidizer_iniital_pressure: float
    fuel_initial_pressure: float
    fuel_pump_pressure_factor: float
    oxidizer_pump_pressure_factor: float
    fuel_pump_specific_power: float
    oxidizer_pump_specific_power: float
    pressurant_gas_constant: float
    pressurant_initial_temperature: float
    oxidizer_pump_efficiency: float
    fuel_pump_efficiency: float
    pressurant_margin_factor: float
    pressurant_tank_structural_factor: float
    propellant_margin_factor: float
    propellant_tanks_structural_factor: float
    ullage_volume_factor: float
    oxidizer_density: float
    fuel_density: float
    propellant_tanks_material_density: float
    pressurant_tank_material_density: float
    propellant_tanks_yield_strength: float
    pressurant_tank_yield_strength: float

    combustion_chamber_material_density: float
    combustion_chamber_yield_strength: float
    combustion_chamber_safety_factor: float
    convergent_throat_bend_ratio: float
    convergent_chamber_bend_ratio: float
    convergent_half_angle: float
    divergent_throat_half_angle: float
    divergent_exit_half_angle: float
    nozzle_type: str
    maximum_wall_temperature: float
    thrust_chamber_wall_emissivity: float
    hot_gas_emissivity: float
    convective_coefficient_mode: str

    recovery_factor: Optional[float] = None
    chamber_throat_area_ratio: Optional[float] = None


class EngineCycle:
    def __init__(self, thrust: float, burn_time: float, combustion_chamber_pressure: float, oxidizer_name: str,
                 fuel_name: str, is_frozen: bool, exit_pressure: float, max_acceleration: float,
                 heat_ratio_pressurant: float, mass_mixture_ratio: float, pressurant_initial_pressure: float,
                 pressurant_final_pressure: float, oxidizer_initial_pressure: float, fuel_initial_pressure: float,
                 fuel_pump_pressure_factor: float, oxidizer_pump_pressure_factor: float,
                 fuel_pump_specific_power: float, oxidizer_pump_specific_power: float, pressurant_gas_constant: float,
                 pressurant_initial_temperature: float, oxidizer_pump_efficiency: float, fuel_pump_efficiency: float,
                 pressurant_margin_factor: float, pressurant_tank_structural_factor: float,
                 propellant_margin_factor: float, tanks_structural_factor: float, ullage_volume_factor: float,
                 oxidizer_density: float, fuel_density: float, tanks_material_density: float,
                 pressurant_tank_material_density: float, tanks_yield_strength: float,
                 pressurant_tank_yield_strength: float, kwak_fix_cycle_type: str,
                 combustion_chamber_material_density: float,
                 combustion_chamber_yield_strength: float, combustion_chamber_safety_factor: float,
                 convergent_half_angle: float, convergent_throat_bend_ratio: float,
                 convergent_chamber_bend_ratio: float,
                 divergent_throat_half_angle: float, divergent_exit_half_angle: float, nozzle_type: str,
                 maximum_wall_temperature: float, thrust_chamber_wall_emissivity: float, hot_gas_emissivity: float,
                 convective_coefficient_mode: str,
                 chamber_throat_area_ratio: Optional[float] = None, recovery_factor: Optional[float] = None
                 , kwak_fix: bool = False):
        self.f_t = thrust  # N
        self.t_b = burn_time  # s
        self.p_cc = combustion_chamber_pressure  # Pa
        self.ox_name = oxidizer_name
        self.fu_name = fuel_name
        self.frozen = is_frozen
        self.p_ex = exit_pressure  # Pa
        self.a_max = max_acceleration  # g0
        self.y_pr = heat_ratio_pressurant  # -
        self.mmr = mass_mixture_ratio  # -
        self.p_pr0 = pressurant_initial_pressure  # Pa
        self.p_pr1 = pressurant_final_pressure  # Pa
        self.p_o_i = oxidizer_initial_pressure  # Pa
        self.p_f_i = fuel_initial_pressure  # Pa
        self.k_fpp = fuel_pump_pressure_factor  # Pa
        self.k_opp = oxidizer_pump_pressure_factor  # Pa
        self.d_fp = fuel_pump_specific_power  # W
        self.d_op = oxidizer_pump_specific_power  # W
        self.r_pr = pressurant_gas_constant  # J/(kg*K)
        self.temp_pr0 = pressurant_initial_temperature  # K
        self.eta_op = oxidizer_pump_efficiency  # -
        self.eta_fp = fuel_pump_efficiency  # -
        self.k_pr = pressurant_margin_factor  # -
        self.k_pr_t = pressurant_tank_structural_factor  # -
        self.k_p = propellant_margin_factor  # -
        self.k_t = tanks_structural_factor  # -
        self.k_u = ullage_volume_factor  # -
        self.rho_o = oxidizer_density  # kg/m3
        self.rho_f = fuel_density  # kg/m3
        self.rho_t = tanks_material_density  # kg/m3
        self.rho_pr_t = pressurant_tank_material_density  # kg/m3
        self.s_y_t = tanks_yield_strength  # Pa
        self.s_y_pr_t = pressurant_tank_yield_strength  # Pa

        self.combustion_chamber_material_density = combustion_chamber_material_density
        self.combustion_chamber_yield_strength = combustion_chamber_yield_strength
        self.combustion_chamber_safety_factor = combustion_chamber_safety_factor
        self.convergent_half_angle = convergent_half_angle
        self.convergent_throat_bend_ratio = convergent_throat_bend_ratio
        self.convergent_chamber_bend_ratio = convergent_chamber_bend_ratio
        self.chamber_throat_area_ratio = chamber_throat_area_ratio
        self.divergent_throat_half_angle = divergent_throat_half_angle
        self.divergent_exit_half_angle = divergent_exit_half_angle
        self.nozzle_type = nozzle_type
        self.maximum_wall_temperature = maximum_wall_temperature
        self.thrust_chamber_wall_emissivity = thrust_chamber_wall_emissivity
        self.hot_gas_emissivity = hot_gas_emissivity
        self.convective_coefficient_mode = convective_coefficient_mode
        self.recovery_factor = recovery_factor
        assert kwak_fix_cycle_type in ['gg', 'ep']
        self.kwak_fix_cycle_type = kwak_fix_cycle_type
        self.kwak_fix = kwak_fix
        self.cstar, self.cf, self.temps, self.transport, self.mw_gamma, self.eps = self.set_cea()
        self.temp_cc, self.temp_th, self.temp_ex = self.temps
        self.mws, self.gammas = self.mw_gamma
        self.mw_cc, self.mw_th, self.mw_ex = self.mws
        self.y_cc, self.y_th, self.y_ex = self.gammas
        self.cps, self.mus, self.ks, self.prandtls = self.transport
        self.cp_cc, self.cp_th, self.cp_ex = self.cps
        self.mu_cc, self.mu_th, self.mu_ex = self.mus
        self.k_cc, self.k_th, self.k_ex = self.ks
        self.prandtl_cc, self.prandtl_th, self.prandtl_k = self.prandtls

    def set_cea(self):
        return get_cea_values(chamber_pressure=self.p_cc, mixture_ratio=self.mmr, exit_pressure=self.p_ex,
                              fuel_name=self.fu_name, ox_name=self.ox_name, isfrozen=self.frozen)

    def reiterate(self):
        self.cstar, self.cf, self.temps, self.transport = self.set_cea()

    @property
    def propellant_mix_name(self):
        if any(name in self.fu_name.lower() for name in ['h2', 'hydrogen']):
            if any(name in self.fu_name.lower() for name in ['gh2', 'gas']):
                fuel = 'GH2'
            elif any(name in self.fu_name.lower() for name in ['lh2', 'liquid']):
                fuel = 'LH2'
        elif any(name in self.fu_name.lower() for name in ['rp1', 'rocket-propellant1', 'rp-1']):
            fuel = 'RP1'

        if any(name in self.ox_name.lower() for name in ['lo2', 'oxygen', 'lox']):
            oxidizer = 'LOX'
        try:
            return f'{oxidizer}/{fuel}'
        except UnboundLocalError:
            warnings.warn(f'Propellant mix name not defined for fu:{self.fu_name}, ox:{self.ox_name}')
            return None

    @property
    def mass_flow(self):
        return self.f_t / (self.cstar * self.cf)

    @property
    def throat_area(self):
        return self.mass_flow * sqrt(constants.R / self.mw_cc * self.temp_cc) / (get_kerckhove(self.y_cc) * self.p_cc)

    @property
    def exit_area(self):
        return self.throat_area * get_expansion_ratio(self.pressure_ratio, self.y_cc)

    @property
    def pressure_ratio(self):
        return self.p_cc / self.p_ex

    @property
    def thrust(self):
        return self.cf * self.cstar * self.mass_flow

    @property
    def simple_specific_impulse(self):
        return self.thrust / self.mass_flow / g

    @property
    def fuel_flow(self):
        return 1 / (self.mmr + 1) * self.mass_flow

    @property
    def oxidizer_flow(self):
        return self.mmr / (self.mmr + 1) * self.mass_flow

    @property
    def oxidizer(self):
        return Propellant(mass_flow=self.oxidizer_flow, burn_time=self.t_b,
                          density=self.rho_o, type='oxidizer', margin_factor=self.k_p)

    @property
    def fuel(self):
        return Propellant(mass_flow=self.fuel_flow, burn_time=self.t_b,
                          density=self.rho_f, type='fuel', margin_factor=self.k_p)

    @property
    def pressurant(self):
        return PressurantTank(oxidizer_volume=self.oxidizer.volume, fuel_volume=self.fuel.volume,
                              fuel_tank_pressure=self.p_f_i, oxidizer_tank_pressure=self.p_o_i, ullage_factor=self.k_u,
                              margin_factor=self.k_pr, structural_factor=self.k_pr_t, initial_pressure=self.p_pr0,
                              final_pressure=self.p_pr1, heat_ratio_pressurant=self.y_pr,
                              gas_constant_pressurant=self.r_pr, initial_temperature_pressurant=self.temp_pr0,
                              material_density=self.rho_pr_t, yield_strength=self.s_y_pr_t)

    @property
    def oxidizer_tank(self):
        return Tank(max_acceleration=self.a_max, ullage_factor=self.k_u, propellant=self.oxidizer,
                    pressurant_tank=self.pressurant, tank_initial_pressure=self.p_o_i,
                    material_density=self.rho_t, yield_strength=self.s_y_t, safety_factor=self.k_t,
                    kwak_fix_cycle_type=self.kwak_fix_cycle_type, kwak_fix=self.kwak_fix)

    @property
    def fuel_tank(self):
        return Tank(max_acceleration=self.a_max, ullage_factor=self.k_u, propellant=self.fuel,
                    pressurant_tank_volume=self.pressurant.tank_volume, tank_initial_pressure=self.p_f_i,
                    material_density=self.rho_t, yield_strength=self.s_y_t, safety_factor=self.k_t,
                    kwak_fix_cycle_type=self.kwak_fix_cycle_type, kwak_fix=self.kwak_fix)

    @property
    def delta_p_oxidizer_pump(self):
        return self.p_cc * self.k_opp - self.p_o_i

    @property
    def delta_p_fuel_pump(self):
        return self.p_cc * self.k_fpp - self.p_f_i

    @property
    def oxidizer_pump(self):
        return Pump(propellant_density=self.oxidizer.density, mass_flow=self.oxidizer.mass_flow,
                    pressure_change=self.delta_p_oxidizer_pump, efficiency=self.eta_op, specific_power=self.d_op)

    @property
    def fuel_pump(self):
        return Pump(propellant_density=self.fuel.density, mass_flow=self.fuel.mass_flow,
                    pressure_change=self.delta_p_fuel_pump, efficiency=self.eta_fp, specific_power=self.d_fp)

    @property
    def combustion_chamber(self):
        return CombustionChamber(
            throat_area=self.throat_area, combustion_chamber_pressure=self.p_cc,
            convergent_half_angle=self.convergent_half_angle, throat_bend_ratio=self.convergent_throat_bend_ratio,
            chamber_bend_ratio=self.convergent_chamber_bend_ratio,
            propellant_mix=self.propellant_mix_name, area_ratio_cc_throat=self.chamber_throat_area_ratio,
            safety_factor=self.combustion_chamber_safety_factor, yield_strength=self.combustion_chamber_yield_strength,
            material_density=self.combustion_chamber_material_density)

    @property
    def nozzle(self):
        return Nozzle(
            throat_area=self.throat_area, nozzle_type=self.nozzle_type, area_ratio=self.eps,
            throat_half_angle=self.divergent_throat_half_angle, exit_half_angle=self.divergent_exit_half_angle
        )

    @property
    def thrust_chamber(self):
        return ThrustChamber(nozzle=self.nozzle, combustion_chamber=self.combustion_chamber,
                             heat_capacity_ratio=self.y_cc)

    @property
    def heat_exchanger(self):
        return HeatExchanger(combustion_temperature=self.temp_cc, combustion_chamber_pressure=self.p_cc,
                             dynamic_viscosity=self.mu_cc, specific_heat_capacity=self.cp_cc, mass_flow=self.mass_flow,
                             maximum_wall_temperature=self.maximum_wall_temperature,
                             thrust_chamber_wall_emissivity=self.thrust_chamber_wall_emissivity,
                             hot_gas_emissivity=self.hot_gas_emissivity, heat_capacity_ratio=self.y_cc,
                             convective_coefficient_mode=self.convective_coefficient_mode,
                             thrust_chamber=self.thrust_chamber, recovery_factor=self.recovery_factor)

    @property
    def pump_power_required(self):
        return self.fuel_pump.power_required + self.oxidizer_pump.power_required

    @property
    def pumps_mass(self):
        return self.fuel_pump.mass + self.oxidizer_pump.mass

    @property
    def tanks_mass(self):
        return self.fuel_tank.mass + self.oxidizer_tank.mass + self.pressurant.tank_mass

    @property
    def props_mass(self):
        return self.oxidizer.mass + self.fuel.mass

    @property
    def mass(self):
        return self.props_mass + self.tanks_mass + self.pumps_mass

    @property
    def mass_ratio(self):
        return (self.mass - self.props_mass) / self.mass

    @property
    def ideal_delta_v(self):
        return self.simple_specific_impulse * log(1 / self.mass_ratio) * g

    @property
    def gravity_delta_v(self, vertical_fraction: float = 0.2):
        return self.ideal_delta_v - constants.g * self.t_b * vertical_fraction

    @property
    def payload_delta_v(self):
        payload = 10  # kg
        return self.simple_specific_impulse * log(self.mass + payload / (self.mass - self.props_mass + payload)) * g

    @property
    def payload_mass_ratio(self):
        if self.mass_u is None:
            ValueError('Payload mass not given, impossible to calculate mass ratio with payload mass')
        return (self.mass + self.mass_u - self.props_mass) / (self.mass + self.mass_u)


if __name__ == '__main__':
    from arguments import base_arguments_o

    throat_area = pi * .05 ** 2
    chamber_pressure = 55e5
    area_ratio = 22
    del base_arguments_o['exit_pressure']
    b = EngineCycle(thrust=100e3,
                    burn_time=160,
                    combustion_chamber_pressure=chamber_pressure,
                    mass_mixture_ratio=5.6,
                    is_frozen=False,
                    fuel_pump_specific_power=15E3,
                    oxidizer_pump_specific_power=20E3,
                    kwak_fix_cycle_type='ep',
                    **base_arguments_o,
                    exit_pressure=.5e5,
                    combustion_chamber_material_density=1,
                    combustion_chamber_yield_strength=1,
                    combustion_chamber_safety_factor=1,
                    convergent_half_angle=radians(30),
                    convergent_throat_bend_ratio=0.8,
                    convergent_chamber_bend_ratio=1.0,
                    chamber_throat_area_ratio=(80 / 50) ** 2,
                    divergent_throat_half_angle=radians(35),
                    divergent_exit_half_angle=radians(5),
                    nozzle_type='bell',
                    maximum_wall_temperature=850,
                    thrust_chamber_wall_emissivity=.8,
                    hot_gas_emissivity=.1,
                    convective_coefficient_mode='Modified Bartz'
                    )

    b.thrust_chamber.show_contour()
    b.thrust_chamber.show_mach()
    b.heat_exchanger.show_adiabatic_wall_temp()
    b.heat_exchanger.show_heat_flux()
    # print(b.throat_area)
    #
    # a = Nozzle(throat_area=throat_area,
    #            nozzle_type='bell',
    #            area_ratio=area_ratio,
    #            throat_half_angle=radians(40),
    #            exit_half_angle=radians(5),
    #            angles_in_rad=True
    #            )
    # d = CombustionChamber(throat_area=throat_area,
    #                       combustion_chamber_pressure=chamber_pressure,
    #                       convergent_half_angle=radians(20),
    #                       characteristic_length=None,
    #                       area_ratio_cc_throat=(80 / 50) ** 2,
    #                       propellant_mix='LOX/LH2',
    #                       chamber_bend_ratio=1.0,
    #                       throat_bend_ratio=.8,
    #                       safety_factor=1,
    #                       material_density=1,
    #                       yield_strength=1)
    # c = ThrustChamber(a, d, 1.207)
    # print(a.length)
    # print(a.radius(a.length), a.exit_radius)
    # print(a.radius(a.length_p), a.radius_p)
    # print(a.radius(0), a.throat_radius)
    # print(b.transport)
    # print(d.length_p)
    # c.show_contour()
    # c.show_mach()
    # print(c.get_mach(0))
    # e = HeatExchanger(
    #     combustion_temperature=3300,
    #     combustion_chamber_pressure=b.p_cc,
    #     dynamic_viscosity=b.mu_cc,
    #     specific_heat_capacity=b.cp_cc,
    #     mass_flow=b.mass_flow,
    #     wall_temperature=850,
    #     chamber_wall_emissivity=.8,
    #     combustion_gas_emissivity=.1,
    #     heat_capacity_ratio=b.y_cc,
    #     convective_coefficient_mode="Modified Bartz",
    #     thrust_chamber=c
    # )
    # print(e.turbulent_recovery_factor, e.tc)
    # e.show_heat_flux_coefficient()
    # e.show_adiabatic_wall_temp()
    # e.show_heat_flux()
    # print(c.surface)
    # print(f'Heat transfer: {e.total_convective_heat_transfer:.3E}')
