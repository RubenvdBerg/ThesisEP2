from EngineCycles.Abstract.EngineCycle import EngineCycle
import pandas as pd
import os


def prettify_attribute(attribute: str):
    words = (word if word != 'secondary' else '2nd' for word in attribute.split('_'))
    return ' '.join(word.capitalize() for word in words)


def get_mass_values(engine: EngineCycle):
    mass_attributes = (
        'fuel_pump',
        'oxidizer_pump',
        'secondary_fuel_pump',
        'turbine',
        'electric_motor',
        'inverter',
        'battery',
        'battery_cooler',
        'injector',
        'thrust_chamber',
        'secondary_exhaust',
        'gas_generator',
        'fuel',
        'oxidizer',
        'pressurant',
        'fuel_tank',
        'oxidizer_tank',
        'pressurant_tank',
        'heat_transfer_section',
    )
    for mass_attribute in mass_attributes:
        component = getattr(engine, mass_attribute, None)
        yield prettify_attribute(mass_attribute), getattr(component, 'mass', None)


def get_mass_distribution(engine: EngineCycle):
    return list(get_mass_values(engine))


def get_mass_distribution_excel(engine: EngineCycle, path='mass_distribution'):
    df = pd.DataFrame(data=get_mass_distribution(engine), columns=['Component', 'Mass [kg]'])
    df.set_index(['Component'])
    with pd.ExcelWriter(path + '.xlsx', engine='xlsxwriter') as writer:
        df.to_excel(writer, header=True, index=False)
    os.system(f'start EXCEL.EXE {path}')


if __name__ == '__main__':
    from EngineArguments import arguments as args
    from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
    from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
    from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
    from EngineCycles.CoolantBleedCycle import CoolantBleedCycle

    design_args = {'thrust': 100e3,
                   'burn_time': 390,
                   'combustion_chamber_pressure': 10e6,
                   'is_frozen': True,
                   'ambient_pressure': None,
                   'exit_pressure_forced': .002e6, }

    cycle_list = (
        (ElectricPumpCycle, args.ep_arguments),
        (GasGeneratorCycle, args.gg_arguments),
        (CoolantBleedCycle, args.cb_arguments),
        (OpenExpanderCycle, args.oe_arguments),
    )

    for Cycle, extra_args in cycle_list:
        complete_args = args.base_arguments | extra_args | design_args | {'verbose': True}
        engine = Cycle(**complete_args)
        get_mass_distribution_excel(engine)
        input('Next?')
