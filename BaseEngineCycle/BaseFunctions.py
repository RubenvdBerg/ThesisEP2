import warnings


def get_propellant_mix_name(fuel_name: str, oxidizer_name: str) -> str:
    if any(name in fuel_name.lower() for name in ['h2', 'hydrogen']):
        if any(name in fuel_name.lower() for name in ['gh2', 'gas']):
            fuel = 'GH2'
        elif any(name in fuel_name.lower() for name in ['lh2', 'liquid']):
            fuel = 'LH2'
        else:
            fuel = 'H2'
    elif any(name in fuel_name.lower() for name in ['rp1', 'rocket-propellant1', 'rp-1']):
        fuel = 'RP1'
    else:
        fuel = None

    if any(name in oxidizer_name.lower() for name in ['lo2', 'oxygen', 'lox']):
        oxidizer = 'LOX'
    else:
        oxidizer = None

    if None in [oxidizer, fuel]:
        warnings.warn(f'Propellant mix name not defined for fu:{fuel_name}, ox:{oxidizer_name}')
    return f'{oxidizer}/{fuel}'


def set_if_not_none(variable, value):
    return value if variable is not None else variable


def only_one_none(a, b, c):
    a, b, c = a is None, b is None, c is None
    return a ^ b ^ c ^ all((a, b, c))

