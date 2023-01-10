propellant_mix_error_message = 'The only accepted propellant combinations are: ["LH2/LOX","LCH4/LOX","RP1/LOX"]'


def get_propellant_mixture(fuel_name: str, oxidizer_name: str) -> str:
    if 'LOX' in oxidizer_name or 'LO2' in oxidizer_name:
        if 'LH2' in fuel_name:
            return 'LH2/LOX'
        elif 'CH4' in fuel_name:
            return 'LCH4/LOX'
        elif 'RP1' in fuel_name:
            return 'RP1/LOX'
    else:
        raise NotImplementedError(propellant_mix_error_message)


def get_characteristic_length(propellant_mix: str) -> float:
    """Find the characteristic length in [m] based on propellant combination."""
    switcher = {'LH2/LOX': 0.89, 'LCH4/LOX': 1.45, 'RP1/LOX': 1.145}
    # R.W. Humble, 1995 - "Space Propulsion Analysis and Design" p.220
    try:
        return switcher[propellant_mix]
    except KeyError:
        raise NotImplementedError(propellant_mix_error_message)


def get_specific_impulse_quality_factor(propellant_mix: str) -> float:
    """Find the specific impulse correction/quality factor [-] based on propellant combination."""
    switcher = {'LH2/LOX': 0.98, 'LCH4/LOX': 0.97, 'RP1/LOX': 0.95}
    try:
        return switcher[propellant_mix]
    except KeyError:
        raise NotImplementedError(propellant_mix_error_message)


def get_mass_mixture_ratio(propellant_mix: str) -> float:
    """Get default mixture ratio [-] based on propellant combination"""
    switcher = {'LH2/LOX': 5.6, 'LCH4/LOX': 3.6, 'RP1/LOX': 2.45}
    try:
        return switcher[propellant_mix]
    except KeyError:
        raise NotImplementedError(propellant_mix_error_message)


def get_initial_propellant_temperature(propellant_name: str) -> float:
    """Get the initial temperature of a propellant in [K]."""
    if 'RP' in propellant_name:
        return 263.6
    elif 'LH2' in propellant_name:
        return 20.25
    elif 'CH4' in propellant_name:
        return 111.0
    elif 'LOX' in propellant_name or 'LO2' in propellant_name:
        return 90.19


def get_prandtl_number_estimate(heat_capacity_ratio: float) -> float:
    return 4 * heat_capacity_ratio / (9 * heat_capacity_ratio - 5)


def get_turbulent_recovery_factor(prandtl_number: float) -> float:
    return prandtl_number ** (1 / 3)
