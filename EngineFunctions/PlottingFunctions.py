from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from EngineFunctions.BaseFunctions import format_attr_name, get_unit

si_prefixes = ('Y', 'Z', 'E', 'P', 'T', 'G', 'M', 'k', '', 'm', '\u03BC', 'n', 'p', 'f', 'a', 'z', 'y')
si_prefix_powers = {key: (i - 8) * 3 for i, key in enumerate(si_prefixes)} | {'h': -2, 'da': -1, 'd': 1, 'c': 2}


def get_si_prefix_power(si_prefix: str):
    try:
        return si_prefix_powers[si_prefix]
    except KeyError:
        raise KeyError(f'Use a valid SI-prefix: {list(si_prefix_powers.keys())}')


def make_axis_string(attribute_name: str, si_prefix: str):
    format_switcher = {
        'power_ratio': r'$\frac{\mathrm{Power\ Mass}}{\mathrm{Burn\ Time}}$',
        'anti_power_ratio': r'$\frac{\mathrm{Anti\ Power\ Mass}}{\mathrm{Burn\ Time}}$',
        'initial_mass_ratio': r'$\frac{\mathrm{Initial\ Mass}}{\mathrm{Burn\ Time}}$',
        'turbine.inlet_flow_state.specific_heat_capacity': r'Turbine Specific Heat Capacity',
        'turbine.inlet_flow_state.heat_capacity_ratio': r'Turbine Heat Capacity Ratio',
        'energy_source_ratio': r'$\frac{\mathrm{Energy\ Source\ Mass}}{\mathrm{Burn\ Time}}$',
        'cc_prop_group_ratio': r'$\frac{\mathrm{CC\ Prop.\ Group\ Mass}}{\mathrm{Burn\ Time}}$',
        'tanks_plus_propellant': r'Tanks + Propellant Mass',
    }
    name = format_switcher[attribute_name] if attribute_name in format_switcher else format_attr_name(attribute_name)
    unit = get_unit(attribute_name)
    return f'{name} [{si_prefix}{unit}]'


def get_class_acronym(EngineClass: EngineCycle):
    if issubclass(EngineClass, ElectricPumpCycle):
        return 'EP'
    elif issubclass(EngineClass, GasGeneratorCycle):
        return 'GG'
    elif issubclass(EngineClass, OpenExpanderCycle):
        return 'OE'


def get_class_color_marker(EngineClass: EngineCycle):
    color_dict = {
        GasGeneratorCycle: ('blue', '^'),
        ElectricPumpCycle: ('green', 's'),
        OpenExpanderCycle: ('red', 'o'),
    }
    return color_dict[EngineClass]

def get_class_from_name(engine_class_name: str) -> EngineCycle:
    cycles = [ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle]
    class_dict = {EngineClass.__name__: EngineClass for EngineClass in cycles}
    return class_dict[engine_class_name]

def adjust_values_to_prefix(values: list, si_prefix: str):
    power = get_si_prefix_power(si_prefix)
    return [val * 10 ** power for val in values]


def adjust_joule_to_watt_hour(values: list):
    return [val /3600 for val in values]


def format_attr_name_for_legend(attribute: str):
    switcher = {
        'turbine.inlet_flow_state.specific_heat_capacity': r'Turbine $c_p$',
        'turbine.inlet_flow_state.heat_capacity_ratio': r'Turbine $\gamma$',
        'tanks_plus_pressurant': r'Tanks + Pressurant',
    }
    if attribute in switcher:
        return switcher[attribute]
    else:
        return format_attr_name(attr_name=attribute)


def format_attr_name_for_axis_label(attribute: str):
    switcher = {
        'turbine.inlet_flow_state.specific_heat_capacity': r'Turbine Specific Heat Capacity',
        'turbine.inlet_flow_state.heat_capacity_ratio': r'Turbine Heat Capacity Ratio',
    }
    if attribute in switcher:
        return switcher[attribute]
    else:
        return format_attr_name(attr_name=attribute)
