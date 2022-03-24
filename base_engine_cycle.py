from math import pi, log
from scipy.constants import g
from cea import get_cstar_cf
from typing import Optional
from scipy.interpolate import interp1d


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


class Injector:
    pass


class CombustionChamber(Structure):
    def __init__(self, propellant_mix: str, throat_area: float, combustion_chamber_pressure,
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
        self.p_cc = combustion_chamber_pressure  # Pa
        super().__init__(**kwargs)

    @property
    def volume(self):
        return self.l_star * self.a_t

    @property
    def mass(self):
        return self.md * self.sf / self.sy * 2 * self.volume * self.p_cc


class ComplexCombustionChamber(CombustionChamber):
    def __init__(self, average_wall_temperature: Optional[float], **kwargs):
        self.t_w_avg = 700 if average_wall_temperature is None else average_wall_temperature  # Kelvin
        super().__init__(**kwargs)
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

    @property
    def volume(self):
        return self.l_star * self.a_t

    @property
    def mass(self):
        return self.md * self.sf / self.sy * 2 * self.volume * self.p_cc


class Nozzle:
    pass


class ThrustChamber:
    pass


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
        self.cstar, self.cf = self.set_cea()

    def set_cea(self):
        return get_cstar_cf(chamber_pressure=self.p_cc, mixture_ratio=self.mmr, exit_pressure=self.p_ex,
                            fuel_name=self.fu_name, ox_name=self.ox_name, isfrozen=self.frozen)

    def reiterate(self):
        self.cstar, self.cf = self.set_cea()

    @property
    def mass_flow(self):
        return self.f_t / (self.cstar * self.cf)

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
    def payload_delta_v(self):
        return self.simple_specific_impulse * log() * g

    @property
    def payload_mass_ratio(self):
        if self.mass_u is None:
            ValueError('Payload mass not given, impossible to calculate mass ratio with payload mass')
        return (self.mass + self.mass_u - self.props_mass) / (self.mass + self.mass_u)
