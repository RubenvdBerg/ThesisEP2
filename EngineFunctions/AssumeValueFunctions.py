def get_characteristic_length(fuel_name: str, oxidizer_name: str) -> float:
    """Find the characteristic length in [m] based on propellant combination."""
    if 'LOX' in oxidizer_name or 'LO2' in oxidizer_name:
        if 'LH2' in fuel_name:
            return 0.89
        elif 'LCH4' in fuel_name:
            return 1.45
        elif 'RP1' in fuel_name:
            return 1.145
    else:
        raise NotImplementedError('The characteristic length can only be estimated for the following propellant '
                                  'combinations: ["LOX/LH2","LOX/LCH4","LOX/RP1"]')


def get_initial_propellant_temperature(propellant_name: str) -> float:
    """Get the initial temperature of a propellant in [K]."""
    if 'RP' in propellant_name:
        return 263.6
    elif 'LH2' in propellant_name:
        return 20.25
    elif 'CH4' in propellant_name:
        return 111.0


def get_prandtl_number_estimate(heat_capacity_ratio:float) -> float:
    return 4 * heat_capacity_ratio / (9 * heat_capacity_ratio - 5)


def get_turbulent_recovery_factor(prandtl_number: float) -> float:
    return prandtl_number ** (1 / 3)

