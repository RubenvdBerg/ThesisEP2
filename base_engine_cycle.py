from math import pi, log, tan, radians, sqrt, cos, sin, asin
from scipy.constants import g, sigma, R
from numpy import roots
from cea import get_cea_values, get_cstar_cf
from typing import Optional
from scipy.interpolate import interp1d
from irt import get_kerckhove, get_expansion_ratio
from dataclasses import dataclass


class Structure:
    def __init__(self, material_density, safety_factor, sigma_yield):
        self.md = material_density  # kg/m3
        self.sf = safety_factor  # -
        self.sy = sigma_yield  # Pa


class Tank(Structure):
    def __init__(self, max_acceleration, ullage_factor, propellant,
                 tank_initial_pressure, kwak_fix_cycle_type, pressurant_tank_volume=None, kwak_fix=False, **kwargs):
        self.ma = max_acceleration  # m/s2 in g0
        self.ku = ullage_factor  # -
        self.pv = propellant.volume  # m3
        self.pd = propellant.density  # kg/m3
        self.pt = propellant.type  # oxidizer or fuel
        self.tip = tank_initial_pressure  # Pa [fuel=.25E6, ox=4E6]
        self.ptv = pressurant_tank_volume  # m3
        self.kwak_fix_cycle_type = kwak_fix_cycle_type
        self.kwak_fix = kwak_fix
        super().__init__(**kwargs)

    @property
    def radius(self):
        return (3 * self.volume / (4 * pi)) ** (1 / 3)

    @property
    def initial_head(self):
        unused_volume = (self.ku - 1) * self.pv
        rest_height = (unused_volume * 3 / (2 * pi)) ** (1 / 3)
        if self.kwak_fix:
            if self.kwak_fix_cycle_type == 'ep':
                if self.pt == 'oxidizer':
                    return 1.91839449096392
                elif self.pt == 'fuel':
                    return 1.65560478870526
            elif self.kwak_fix_cycle_type == 'gg':
                if self.pt == 'oxidizer':
                    return 1.92539045861846
                elif self.pt == 'fuel':
                    return 1.71033897378923
        return 2 * self.radius - rest_height

    @property
    def total_lower_pressure(self):
        return self.tip + self.pd * self.ma * g * self.initial_head

    @property
    def total_upper_pressure(self):
        return self.tip + self.pd * self.ma * g * (self.initial_head - self.radius)

    @property
    def mass(self):  # Spherical Tanks
        return self.sf * 3 * self.md / (4 * self.sy) * self.volume * (
                self.total_upper_pressure + self.total_lower_pressure)

    @property
    def volume(self):
        if self.pt == 'oxidizer':
            return self.pv * self.ku + self.ptv
        elif self.pt == 'fuel':
            return self.pv * self.ku


class Propellant:
    def __init__(self, mass_flow, burn_time, density, propellant_type, margin_factor=1.01):
        # LOX density = 1126.1
        # RP-1 density = 804.2
        types = {'oxidizer', 'fuel'}
        if propellant_type not in types:
            raise ValueError(f'propellant_type must be one of {types}')
        self.type = propellant_type  # oxidizer or fuel
        self.mass_flow = mass_flow  # kg/s
        self.burn_time = burn_time  # s
        self.density = density  # kg/m3
        self.margin_factor = margin_factor  # -

    @property
    def mass(self):
        return self.mass_flow * self.burn_time * self.margin_factor

    @property
    def volume(self):
        return self.mass / self.density

    @property
    def volumetric_flow(self):
        return self.mass_flow / self.density


class PressurantTank(Structure):
    def __init__(self, oxidizer_volume, fuel_volume, fuel_tank_pressure, oxidizer_tank_pressure,
                 ullage_factor, margin_factor, structural_factor, initial_pressure,
                 final_pressure, heat_ratio_pressurant, gas_constant_pressurant, initial_temperature_pressurant,
                 **kwargs):
        self.ov = oxidizer_volume  # m3
        self.fv = fuel_volume  # m3
        self.ftp = fuel_tank_pressure  # Pa
        self.otp = oxidizer_tank_pressure  # Pa
        self.ku = ullage_factor  # -
        self.km = margin_factor  # -
        self.kt = structural_factor  # -
        self.p0 = initial_pressure  # Pa
        self.p1 = final_pressure  # Pa
        self.hy = heat_ratio_pressurant  # - [specific heat ratio of helium (helium gamma -> hy)]
        self.Rh = gas_constant_pressurant  # J/(kg*K)
        self.Th0 = initial_temperature_pressurant  # K
        super().__init__(safety_factor=structural_factor, **kwargs)

    @property
    def mass(self):
        return (self.km * self.ku * self.hy / (self.Rh * self.Th0) * (self.otp * self.ov + self.ftp * self.fv)
                / (1 - (self.p1 / self.p0)))

    @property
    def tank_volume(self):
        return self.mass * self.Rh * self.Th0 / self.p0

    @property
    def tank_mass(self):
        return self.kt * 3 / 2 * self.tank_volume * self.p0 * self.md / self.sy


class Pump:
    def __init__(self, propellant_density: float, mass_flow: float, pressure_change: float, efficiency: float,
                 specific_power: float):
        self.rho_p = propellant_density  # kg/m3
        self._mass_flow = mass_flow  # kg/s
        self.delta_p = pressure_change  # Pa
        self.eta = efficiency  # -
        self.sp = specific_power  # W/kg

    @property
    def mass_flow(self):
        return self._mass_flow

    @mass_flow.setter
    def mass_flow(self, x):
        self._mass_flow = x

    @property
    def q_p(self):
        return self.mass_flow / self.rho_p

    @property
    def power_required(self):
        return self.q_p * self.delta_p / self.eta

    @property
    def mass(self):
        return self.power_required / self.sp


class Injector(Structure):
    def __init__(self, combustion_chamber_area, combustion_chamber_pressure, propellant_is_gas=False, **kwargs):
        self.a_cc = combustion_chamber_area
        self.p_cc = combustion_chamber_pressure
        self.is_gas = propellant_is_gas
        super().__init__(**kwargs)

    @property
    def pressure_drop(self):
        # Mota 2008 -> Kesaev and Almeida 2005
        f = .4 if self.is_gas else .8
        return f * 10E2 * sqrt(10 * self.p_cc)

    @property
    def thickness(self):
        # Zandbergen -> 4 times pressure drop, ASME code for welded flat caps on internal pressure vessels
        diameter = 2 * sqrt(self.a_cc / pi)
        pressure = 4 * self.pressure_drop
        c = 6 * 3.3 / 64
        return diameter * sqrt(c * pressure / self.sy) * self.sf

    @property
    def mass(self):
        return self.a_c * self.thickness * self.md


class CombustionChamber(Structure):
    def __init__(self, throat_area: float, combustion_chamber_pressure, convergent_half_angle: float,
                 in_radians: bool = False, throat_bend_ratio: float = 1., chamber_bend_ratio: float = .8,
                 propellant_mix: Optional[str] = None,
                 characteristic_length: Optional[float] = None, **kwargs):

        # Set default value for characteristic length
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
                    f'The specified propellant mix {propellant_mix} does not have a default characteristic length for the combustion chamber, specify one manually or select a propellant mix that does [{options.keys}]')
        else:
            self.l_star = characteristic_length
        # Initialize other parameters
        self.a_t = throat_area  # m2
        self.r_t = sqrt(self.a_t / pi)  # m
        self.p_cc = combustion_chamber_pressure  # Pa
        self.cc_bend = chamber_bend_ratio
        self.th_bend = throat_bend_ratio
        self.ha_conv = convergent_half_angle if in_radians else radians(convergent_half_angle)
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
        # Humble 1995 p.222
        nozzle_throat_diameter = 2 * sqrt(self.a_t / pi)
        ac_at = 8.0 * nozzle_throat_diameter ** 2.4 + 1.25
        return ac_at * self.a_t

    @property
    def radius(self):
        return sqrt(self.area / pi)

    @property
    def convergent_length(self):
        return self.length_p + self.r_a * sin(self.ha_conv)

    @property
    def length_p(self):
        z = self.radius - self.r_t - (self.r_u + self.r_a) * (1 - cos(self.ha_conv))
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


class Nozzle:
    def __init__(self, throat_area, nozzle_type, area_ratio, throat_divergence_half_angle, exit_divergence_half_angle,
                 angles_in_rad=True):
        self.a_t = throat_area
        self.eps = area_ratio
        self.type = nozzle_type
        self.half_angle_th = throat_divergence_half_angle if angles_in_rad else radians(throat_divergence_half_angle)
        self.half_angle_ex = exit_divergence_half_angle if angles_in_rad else radians(exit_divergence_half_angle)
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
            a = 10
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
            radius = roots([self.a, self.b, self.c - distance_from_throat])[0]
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


class ThrustChamber:
    def __init__(self, nozzle: type(Nozzle), combustion_chamber: type(CombustionChamber)):
        self.nozzle = nozzle
        self.cc = combustion_chamber

    def radius(self, distance_from_throat):
        if distance_from_throat < -self.cc.total_length:
            raise ValueError(
                f'Radius of thrust chamber cannot be calculated before the injector plate ({self.cc.total_length:.4e} m before the throat)')
        elif distance_from_throat < -self.cc.convergent_length:
            return self.cc.radius
        elif distance_from_throat < 0:
            return self.cc.convergent_radius(-distance_from_throat)
        else:
            return self.nozzle.radius(distance_from_throat)

    def show_contour(self):
        distances = linspace(-self.cc.convergent_length, self.nozzle.length, 1000)
        plt.plot(distances, [self.radius(x) for x in distances])
        plt.show()


class HeatExchanger:
    def __init__(self, combustion_temperature, combustion_chamber_pressure, transport_properties, contour_function,
                 mass_flow):
        self.pcc = combustion_chamber_pressure
        self.tc = combustion_temperature
        self.eps_cw = emissivity_chamber_wall
        self.eps_cg = emissivity_combustion_gas
        self.m_flow = mass_flow
        self.tw = wall_temperature
        self.f_contour = contour_radius_function
        self.specific_heat_capacities, _, self.viscosities, self.prandtl_numbers = transport_properties

    @property
    def mass(self):
        return

    @property
    def pressure_loss(self):
        return .15 * self.pcc

    @property
    def temp_diff(self):
        return

    @property
    def netto_average_wall_radiative_heat_flux(self):  # q_rad [W/m2]
        # Heat Transfer Handbook, A. Benjan 2003, Eq. 8.69
        return sigma(self.tc ** 4 - self.tw ** 4) / (1 / self.eps_cg + (1 / self.eps_cw) - 1)

    @property
    def prandtl(self):
        return 4 * self.y / (9 * self.y - 5)

    def get_convective_heat_flux(self, distance_from_throat):
        diameter = 2 * self.f_contour(distance_from_throat)
        return convective_heat_flux(
            mode="Cornellisse",
            mass_flow=self.m_flow,
            diameter=diameter,
            dynamic_viscosity=self.viscosities[1],
            specific_heat_capacity_p=self.specific_heat_capacities[1],
            prandtl_number=self.prandtl_numbers[1],
            )

    # @property
    # def convective_heat_flux(self):
    #     if self.convective_mode == "Modified Bartz":
    #         return 0.026 * 1.213 * self.mass_flow**.8 * self.diameter**-1.8 * self.mu**.2 * self.cp * self.prandtl**-.6 * (self.tc / self.tf)**.68
    #     elif self.convective_mode == "Cornellisse":
    #         return 0.023 * 1.213 * self.mass_flow**.8 * self.diameter**-1.8 * self.mu**.2 * self.cp * self.prandtl**-2/3
    #     elif self.convective_mode == "Standard Bartz":
    #         raise NotImplementedError("Convective heat transfer for the standard bartz equation has not been implemented")
    #     else:
    #         raise ValueError("Improper convective_mode given for calculation of the convective heat transfer")


def convective_heat_flux(mode: str, mass_flow: float, diameter: float, dynamic_viscosity: float,
                         specific_heat_capacity_p: float, prandtl_number: float, total_temperature: float = None,
                         film_temperature: float = None):
    # Zandbergen 2017 p.161
    if mode == "Modified Bartz":
        return 0.026 * 1.213 * mass_flow ** .8 * diameter ** -1.8 * dynamic_viscosity ** .2 * specific_heat_capacity_p * prandtl_number ** -.6 * (
                total_temperature / film_temperature) ** .68
    elif mode == "Cornellisse":
        return 0.023 * 1.213 * mass_flow ** .8 * diameter ** -1.8 * dynamic_viscosity ** .2 * specific_heat_capacity_p * prandtl_number ** -2 / 3
    elif mode == "Standard Bartz":
        raise NotImplementedError("Convective heat transfer for the standard bartz equation has not been implemented")
    else:
        raise ValueError(
            "Improper convective_mode given for calculation of the convective heat transfer, pick one of [Modified Bartz, Cornellise, Standard Bartz]")


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
                 pressurant_tank_yield_strength: float, kwak_fix_cycle_type: str, kwak_fix: bool = False):
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
        assert kwak_fix_cycle_type in ['gg', 'ep']
        self.kwak_fix_cycle_type = kwak_fix_cycle_type
        self.kwak_fix = kwak_fix
        self.cstar, self.cf, self.temps, self.transport, self.mw_gamma, self.eps = self.set_cea()
        self.temp_cc, self.temp_th, self.temp_ex = self.temps
        self.mws, self.gammas = self.mw_gamma
        self.mw_cc, self.mw_th, self.mw_ex = self.mws
        self.y_cc, self.y_th, self.y_ex = self.gammas

    def set_cea(self):
        return get_cea_values(chamber_pressure=self.p_cc, mixture_ratio=self.mmr, exit_pressure=self.p_ex,
                              fuel_name=self.fu_name, ox_name=self.ox_name, isfrozen=self.frozen)

    def reiterate(self):
        self.cstar, self.cf, self.temps, self.transport = self.set_cea()

    @property
    def mass_flow(self):
        return self.f_t / (self.cstar * self.cf)

    @property
    def throat_area(self):
        return self.mass_flow * sqrt(R / self.mw_cc * self.temp_cc) / (get_kerckhove(self.y_cc) * self.p_cc)

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
                          density=self.rho_o, propellant_type='oxidizer', margin_factor=self.k_p)

    @property
    def fuel(self):
        return Propellant(mass_flow=self.fuel_flow, burn_time=self.t_b,
                          density=self.rho_f, propellant_type='fuel', margin_factor=self.k_p)

    @property
    def pressurant(self):
        return PressurantTank(oxidizer_volume=self.oxidizer.volume, fuel_volume=self.fuel.volume,
                              fuel_tank_pressure=self.p_f_i, oxidizer_tank_pressure=self.p_o_i, ullage_factor=self.k_u,
                              margin_factor=self.k_pr, structural_factor=self.k_pr_t, initial_pressure=self.p_pr0,
                              final_pressure=self.p_pr1, heat_ratio_pressurant=self.y_pr,
                              gas_constant_pressurant=self.r_pr, initial_temperature_pressurant=self.temp_pr0,
                              material_density=self.rho_pr_t, sigma_yield=self.s_y_pr_t)

    @property
    def oxidizer_tank(self):
        return Tank(max_acceleration=self.a_max, ullage_factor=self.k_u, propellant=self.oxidizer,
                    pressurant_tank_volume=self.pressurant.tank_volume, tank_initial_pressure=self.p_o_i,
                    material_density=self.rho_t, sigma_yield=self.s_y_t, safety_factor=self.k_t,
                    kwak_fix_cycle_type=self.kwak_fix_cycle_type, kwak_fix=self.kwak_fix)

    @property
    def fuel_tank(self):
        return Tank(max_acceleration=self.a_max, ullage_factor=self.k_u, propellant=self.fuel,
                    pressurant_tank_volume=self.pressurant.tank_volume, tank_initial_pressure=self.p_f_i,
                    material_density=self.rho_t, sigma_yield=self.s_y_t, safety_factor=self.k_t,
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
        return self.ideal_delta_v - g * self.t_b * vertical_fraction

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
    import matplotlib.pyplot as plt
    from numpy import linspace

    b = EngineCycle(thrust=6770e3,
                    burn_time=160,
                    combustion_chamber_pressure=7E6,
                    mass_mixture_ratio=2.27,
                    is_frozen=False,
                    fuel_pump_specific_power=15E3,
                    oxidizer_pump_specific_power=20E3,
                    kwak_fix_cycle_type='ep',
                    **base_arguments_o)
    # print(b.throat_area)
    a = Nozzle(throat_area=b.throat_area,
               nozzle_type='bell',
               area_ratio=16,
               throat_divergence_half_angle=radians(45),
               exit_divergence_half_angle=radians(15)
               )
    d = CombustionChamber(throat_area=b.throat_area,
                          combustion_chamber_pressure=7E6,
                          convergent_half_angle=30,
                          propellant_mix='LOX/RP1', safety_factor=1, material_density=1, sigma_yield=1)
    c = ThrustChamber(a, d)
    print(sqrt(b.exit_area / pi), sqrt(b.throat_area / pi))
    print(a.length)
    print(a.radius(a.length), a.exit_radius)
    print(a.radius(a.length_p), a.radius_p)
    print(a.radius(0), a.throat_radius)
    c.show_contour()
    e = HeatExchanger()
