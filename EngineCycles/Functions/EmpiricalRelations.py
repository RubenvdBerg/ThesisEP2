from math import log, tanh, sqrt, pi
import warnings
from typing import Optional


def get_friction_factor(roughness_height: float, reynolds_number: float, diameter: float) -> float:
    """Estimate friction factor for a circular pipe.

    Theory can be found in:
    "New correlations of single-phase friction factor for turbulent pipe flow and evaluation of existing single-phase
    friction factor correlations" - Fang et al. 2011
    """
    e = roughness_height
    D = diameter
    re = reynolds_number
    rr = e / D
    fds = {}

    # Ranges overlap such that there is a smooth transition (See also ISBN: 9780081024874)
    if 3000 < re and 0 < rr <= .05:
        if re > 1e8:
            warnings.warn(
                'Explicit Colebrook Correlation for Friction Factor used with Reynolds number higher than verified range (above 1e8)')
        # Explicit Colebrook Correlation by Fang et al 2011
        fd1 = 1.1613 * (log(.234 * rr ** 1.1007 - 60.525 * re ** -1.1105 + 56.291 * re ** 1.0712, 10)) ** -2
        fds['Fang'] = fd1
    if rr > 0:
        # Fall back to Nikuradse correlation for turbulent flow in non-smooth pipes
        fd2 = 8 * (2.457 * log(3.707 * (1 / rr))) ** -2
        fds['Nikuradse'] = fd2
    if 1e7 > re > 1e3:
        # Fall back to Blasius equation for turbulent flow in smooth pipes
        fd3 = .184 * re ** -.2
        fds['BlasiusHigh'] = fd3
    if 1e5 > re > 1000:
        # Blasius again
        fd4 = .316 * re ** -.25
        fds['BlasiusLow'] = fd4
    if 3000 > re > 0:
        # Fall back to Hagen-Pouseille law for internal laminar flow
        fd5 = 64 / re
        fds['Laminar'] = fd5
    return max(fds.values())


def get_coolant_convective_heat_transfer_coeff(coolant_conductivity: float,
                                               characteristic_dimension: float,
                                               coolant_reynolds: float,
                                               coolant_prandtl_number: float,
                                               coolant_bulk_temp: float,
                                               wall_temp: float,
                                               mode: str = 'SiederTate',
                                               length_to_start_channel: Optional[float] = None, ) -> float:
    k_c = coolant_conductivity
    D = characteristic_dimension
    re = coolant_reynolds
    pr = coolant_prandtl_number
    t_b = coolant_bulk_temp
    t_w = wall_temp

    if mode == 'Taylor':
        a = 0.023
        exponent = 0.57 - 1.59 * D / length_to_start_channel
    elif mode == 'SiederTate':
        a = 0.025
        exponent = .55
    elif mode == 'DittusBoelter':
        a = 0.023
        exponent = 0
    else:
        raise ValueError(
            'Please select a valid coolant heat transfer coefficient mode: [Taylor, SiederTate, DittusBoelter]')

    return a * k_c / D * re ** .8 * pr ** .4 * (t_b / t_w) ** exponent


def get_roughness_correction(bulk_prandtl_number: float, bulk_reynolds_number: float, roughness_height: float,
                             diameter: float) -> float:
    """Return factor which corrects the convective heat transfer coefficient to account for flow channel roughness."""
    pr = bulk_prandtl_number
    re = bulk_reynolds_number

    fd = get_friction_factor(roughness_height=roughness_height, reynolds_number=re, diameter=diameter)
    fd0 = get_friction_factor(roughness_height=0, reynolds_number=re, diameter=diameter)
    k = fd / fd0
    b = 1.5 * pr ** (-1. / 6.) * re ** (-1. / 8.)
    c_rough = ((1 + b * (pr - 1)) / (1 + b * (pr * k - 1)) * k)
    return c_rough


def get_fin_correction(convective_heat_transfer_coeff: float, fin_thickness: float, channel_height: float, channel_width: float, wall_conductivity: float):
    """Return factor which corrects the convective heat transfer coefficient to account for conduction through fins in
    between flow channels.
    """
    # Fin Correction: Luka Denies - Regenerative cooling analysis of oxygen/methane rocket engines 2015
    h_c = convective_heat_transfer_coeff
    th_fin = fin_thickness
    ht_c = channel_height
    w_c = channel_width
    k = wall_conductivity
    b = (2 * h_c * th_fin / k) ** .5 / th_fin * ht_c
    eta = tanh(b) / b
    c_fin = (w_c + eta * 2 * ht_c) / (w_c + th_fin)
    return c_fin


def get_hot_gas_convective_heat_transfer_coefficient(
        mass_flow: float, local_diameter: float, dynamic_viscosity: float, specific_heat_capacity: float,
        prandtl_number: float, stagnation_temp: float, film_temp: float, mode: str = 'ModifiedBartz',
        local_mach: Optional[float] = None, wall_temp: Optional[float] = None,
        heat_capacity_ratio: Optional[float] = None, throat_radius_of_curvature: Optional[float] = None,
        combustion_chamber_pressure: Optional[float] = None, characteristic_velocity: Optional[float] = None,
        nozzle_throat_diameter: Optional[float] = None
):
    mf = mass_flow
    di = local_diameter
    mu = dynamic_viscosity
    cp = specific_heat_capacity
    pr = prandtl_number
    t0 = stagnation_temp
    tf = film_temp
    if mode == 'Bartz':
        m = local_mach
        tw = wall_temp
        y = heat_capacity_ratio
        dt = nozzle_throat_diameter
        rc = throat_radius_of_curvature
        p0 = combustion_chamber_pressure
        cstar = characteristic_velocity
        if any([prop is None for prop in [m, tw, y, dt, rc, p0, cstar]]):
            raise ValueError('Not all optional parameters, which are required for [Bartz]-mode, are specified.')
        s = (1 + m ** 2 * (y - 1) / 2) ** -.12 / (.5 + .5 * (tw / t0) * (1 + m ** 2 * (y - 1) / 2)) ** .68
        h_hg = (0.026 * dt ** -0.2 * mu ** 0.2 * cp * pr ** -0.6 * (p0 / cstar) ** 0.8
                * (dt / rc) ** .1 * (dt / di) ** 1.8 * s)
    elif mode == 'ModifiedBartz':
        h_hg = 0.026 * 1.213 * mf ** .8 * di ** -1.8 * mu ** .2 * cp * pr ** -.6 * (t0 / tf) ** .68
    else:
        raise ValueError('Select proper mode for estimation of the hot gas convective heat transfer coefficient')
    return h_hg


def get_film_temp(wall_temp:float, adiabatic_wall_temp:float, static_temp:float) -> float:
    """Return film temperature as defined by Ziebland et al. 1971 - Heat Transfer in Rocket Engines"""
    return static_temp * .28 + wall_temp * .5 + adiabatic_wall_temp * .22


def get_chamber_throat_area_ratio_estimate(throat_area: float):
    # Humble 1995 p.222
    throat_diameter_in_cm = 2 * sqrt(throat_area / pi) * 1e2
    area_ratio_chamber_throat = (8.0 * throat_diameter_in_cm ** -.6 + 1.25)
    return area_ratio_chamber_throat