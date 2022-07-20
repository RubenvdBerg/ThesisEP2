from scipy.constants import g, gas_constant, Boltzmann, pi, R
from typing import Optional
from scipy.optimize import fsolve
from math import sqrt
import numpy as np
import warnings


# from sympy import nsolve, solve, Symbol, sqrt


def get_mass_flow(chamber_pressure: float, throat_area: float, chamber_temperature: float, molecular_weight: float,
                  specific_heat_ratio: float) -> float:
    # Molecular Weight has g/mol for SI units, but mass flow should be kg/s
    r = gas_constant / molecular_weight * 10E3
    return get_kerckhove(specific_heat_ratio) * chamber_pressure * throat_area / sqrt(r * chamber_temperature)


def get_expansion_ratio(mach_number: float, heat_capacity_ratio: float):
    m = mach_number
    y = heat_capacity_ratio
    return ((y + 1) / 2) ** -((y + 1) / (2 * (y - 1))) * ((1 + (y - 1) / 2 * m ** 2) ** ((y + 1) / (2 * (y - 1)))) / m


def get_kerckhove(heat_capacity_ratio: float) -> float:
    y = heat_capacity_ratio
    return np.sqrt(y) * (2 / (y + 1)) ** ((y + 1) / (2 * (y - 1)))


def get_expansion_ratio_from_p_ratio(pressure_ratio: float, specific_heat_ratio: float) -> float:
    y = specific_heat_ratio
    G = get_kerckhove(y)
    pe_pc = pressure_ratio ** -1
    p1 = (2 * y / (y - 1)) * pe_pc ** (2 / y) * (1 - pe_pc ** ((y - 1) / y))
    return float(G / np.sqrt(p1))


def get_pressure_ratio(expansion_ratio: float, specific_heat_ratio: float, sympy_solve: bool = False) -> float:
    pr = Symbol('pr')
    pe_pc = pr ** -1
    y = specific_heat_ratio
    gamma = get_kerckhove(y)
    p1 = (2 * y / (y - 1)) * pe_pc ** (2 / y)
    p2 = (1 - pe_pc ** ((y - 1) / y))
    equation = gamma / np.sqrt(p1 * p2) - expansion_ratio
    if sympy_solve:
        return solve(equation, pr, implicit=True)
    else:
        return nsolve(equation, pr, 100)


def get_pressure_ratio_fsolve(expansion_ratio: float, specific_heat_ratio: float, guess: float = 2.) -> float:
    def func(x):
        eps = float(get_expansion_ratio_from_p_ratio(x, specific_heat_ratio))
        return np.array(eps - expansion_ratio, dtype=float)

    solution = float(fsolve(func, np.array(guess, dtype=float))[0])
    return solution
    # # Check
    # if check:
    #     if not np.isclose(func(solution), [0.0]):
    #         message = f'Found pressure ratio does not match given expansion ratio. ' \
    #                   f'Given:[{expansion_ratio}], Found:[{get_expansion_ratio(solution, specific_heat_ratio)}]'
    #         raise ValueError(message)


def get_characteristic_velocity(molecular_weight: float, chamber_temperature: float,
                                specific_heat_ratio: float) -> float:
    # Molecular Weight has g/mol for SI units, but gas constant in J / (mol * K) -> m^2 * kg / (mol * K * s^2)
    r = gas_constant / molecular_weight * 1E3  # J / (kg * K)
    gamma = get_kerckhove(specific_heat_ratio)
    return 1 / gamma * sqrt(r * chamber_temperature)


def get_thrust_coefficient(pressure_ratio: float, specific_heat_ratio: float, expansion_ratio: float,
                           chamber_pressure: Optional[float] = 1, ambient_pressure: Optional[float] = 0,
                           ideal_expansion: bool = True) -> float:
    pe_pc = 1 / pressure_ratio
    y = specific_heat_ratio
    if ideal_expansion:
        pressure_term = 0
    else:
        pressure_term = (pe_pc - ambient_pressure / chamber_pressure) * expansion_ratio
    return get_kerckhove(y) * sqrt(2 * y / (y - 1) * (1 - pe_pc ** ((y - 1) / y))) + pressure_term


def get_specific_impulse(thrust_coefficient: float, characteristic_velocity: float) -> float:
    return thrust_coefficient * characteristic_velocity / g


def thrust_force(specific_impulse: float, mass_flow: float) -> float:
    return specific_impulse * mass_flow * g


def area_to_radius(area):
    return sqrt(area / pi)


def get_exhaust_velocity(molecular_weight: float, specific_heat_ratio: float, chamber_temperature: float,
                         pressure_ratio: float) -> float:
    y = specific_heat_ratio
    pe_pc = 1 / pressure_ratio
    return sqrt(2 * y / (y - 1) * gas_constant / molecular_weight * chamber_temperature * (1 - pe_pc ** ((y - 1) / y)))


def is_choked(pressure_ratio: float, specific_heat_ratio: float):
    pcrit_pc = (2 / (specific_heat_ratio + 1)) ** (specific_heat_ratio / (specific_heat_ratio - 1))
    pe_pc = 1 / pressure_ratio
    return pe_pc < pcrit_pc


def get_isp_simple_vac(molecular_weight: float, specific_heat_ratio: float, chamber_temperature: float,
                       chamber_pressure: float, expansion_ratio: float, complete: bool = False):
    pr = get_pressure_ratio(expansion_ratio, specific_heat_ratio)
    c_star = get_characteristic_velocity(molecular_weight, chamber_temperature, specific_heat_ratio)
    c_f = get_thrust_coefficient(pr, specific_heat_ratio, expansion_ratio, chamber_pressure, ideal_expansion=False)
    i_sp = get_specific_impulse(c_f, c_star)
    if complete:
        return molecular_weight, specific_heat_ratio, chamber_pressure, c_star, c_f, pr, i_sp
    return i_sp


def get_mp_from_isp_itot(total_impulse: float, specific_impulse: float) -> float:
    # Specific Impulse (I_sp) in s, Total Impulse (I_tot) in Ns, g0 in m/s2
    # Assumption: Thrust is constant
    return total_impulse / specific_impulse / g


def get_knudsen(temperature: float, particle_shell_diameter: float, total_pressure: float,
                length_scale: float) -> float:
    return Boltzmann * temperature / (sqrt(2) * pi * particle_shell_diameter ** 2 * total_pressure * length_scale)


def get_mean_free_path(dynamic_viscosity: float, pressure: float, temperature: float, molar_mass: float) -> float:
    # Note: Molar Mass -> kg/mol
    return dynamic_viscosity / pressure * sqrt(pi * R * temperature / (2 * molar_mass))


def get_reynolds(dynamic_viscosity: float, density: float, flow_velocity: float, characteristic_length: float) -> float:
    return density * flow_velocity * characteristic_length / dynamic_viscosity


def get_sonic_velocity(specific_heat_ratio: float, molar_mass: float, temperature: float):
    return sqrt(specific_heat_ratio * (R / molar_mass) * temperature)


def get_local_mach(local_area_ratio, is_subsonic=False, heat_capacity_ratio=1.14):
    a_at = local_area_ratio
    y = heat_capacity_ratio
    p = 2 / (y + 1)
    q = 1 - p
    if is_subsonic:
        r, a, s = a_at ** 2, p ** (1 / q), 1
    else:
        r, a, s = a_at ** (2 * q / p), q ** (1 / p), -1
    r2 = (r - 1) / (2 * a)
    initial_guess = 1 / ((1 + r2) + sqrt(r2 * (r2 + 2)))

    # Python version B4Wind Method by Dennis Yoder from NASA

    def get_f_and_derivs(x):
        f = (p + q * x) ** (1 / q) - r * x
        df = (p + q * x) ** (1 / q - 1) - r
        ddf = p * ((p + q * x) ** (1 / q - 2))
        return f, df, ddf

    def newton_raphson_plus(x):
        while True:
            f, df, ddf = get_f_and_derivs(x)
            xnew = x - 2 * f / (df - sqrt(df ** 2 - 2 * f * ddf))
            if abs(xnew - x) / xnew < .001:
                break
            x = xnew
        return xnew

    final = newton_raphson_plus(initial_guess)
    return sqrt(final) ** s
