from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle, GasGeneratorCycle_DoubleTurbine, \
    GasGeneratorCycle_DoubleTurbineSeries
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle, OpenExpanderCycle_DoubleTurbine
from EngineCycles.CoolantBleedCycle import CoolantBleedCycle
from EngineComponents.Abstract.FlowState import FlowState, ManualFlowState
from EngineComponents.Base.Pump import Pump
from EngineComponents.Other.Battery import Battery
from EngineComponents.Abstract.ElectricalComponent import ElectricalComponent
from EngineComponents.Other.Turbine import Turbine
from numpy import isclose
from typing import Optional, Iterator
from PIL import Image, ImageDraw, ImageFont
from EngineFunctions.BaseFunctions import format_si




def make_mass_schematic(engine: EngineCycle):
    name_switcher = {
        ElectricPumpCycle: ('EP', (1100, 500), (5, 9, 10, 11, 12,)),
        GasGeneratorCycle: ('GG', (1100, 500), (9, 12, 13, 14, 15, 16)),
        CoolantBleedCycle: ('CB', (1100, 500), (9, 11, 12, 13, 14, 15, 16)),
        OpenExpanderCycle: ('OE', (1250, 500), (9, 11, 13, 14, 15, 16)),
        GasGeneratorCycle_DoubleTurbine: ('GG2', (0, 0), (9, 12, 13, 14, 15, 16)),
        OpenExpanderCycle_DoubleTurbine: ('OE2', (0, 0), (9, 11, 13, 14, 15, 16)),
        GasGeneratorCycle_DoubleTurbineSeries: ('GG3', (0, 0), (9, 12, 13, 14, 15, 16)),
    }
    for base_class, info in name_switcher.items():
        if issubclass(type(engine), base_class):
            cycle_name, start_coord, pop_tuple = info

    values = get_values(engine)
    strings = list(format_values(values))
    strings = [string for (i, string) in enumerate(strings, start=1) if i not in pop_tuple]
    fontsize = 63
    font_file = r'C:\Users\rvand\PycharmProjects\ThesisEP2\plots\Imaging\CamingoCode-Regular.ttf'
    myfont = ImageFont.truetype(font_file, fontsize)
    image_path = rf'C:\Users\rvand\PycharmProjects\ThesisEP2\plots\Imaging\MassSchematics\{cycle_name}_Mass_Clean.png'

    totals = [f'Initial: {engine.initial_mass:>8.1f} kg', f'Final  : {engine.final_mass:>8.1f} kg']
    start_x = start_coord[0]
    start_y = start_coord[1]
    dy = fontsize * 1.3
    with Image.open(image_path) as img:
        drawer = ImageDraw.Draw(img)
        for i, string in enumerate(strings):
            coordinate = [start_x, start_y + i*dy]
            drawer.text(coordinate, string , fill=(0, 0, 0), font=myfont)
        for i, string in enumerate(totals):
            coordinate = [start_x + 450, start_y - 300 + i*dy]
            drawer.text(coordinate, string , fill=(0, 0, 0), font=myfont)
        img.show()


def get_values(engine: EngineCycle):
    mass_attributes = (
        'fuel_tank',
        'oxidizer_tank',
        'fuel_pump',
        'oxidizer_pump',
        'turbine',
        'heat_transfer_section',
        'injector',
        'thrust_chamber',
        'splitter',
        'secondary_exhaust',
        'gas_generator',
        'secondary_fuel_pump',
        'electric_motor',
        'inverter',
        'battery',
        'battery_cooler',
        'pressurant_tank',
        'pressurant',
        'fuel',
        'oxidizer',
    )

    mass_components = [getattr(engine, attribute, None) for attribute in mass_attributes]
    for component in mass_components:
        getattr(component, 'mass', 0)
    return [getattr(component, 'mass', 0) for component in mass_components]


def format_values(values: list[float, ...]) -> Iterator[str]:
    mass_names = (
        'Fuel Tank',
        'Oxidizer Tank',
        'Fuel Pump',
        'Oxidizer Pump',
        'Turbine',
        'Heat Exchanger',
        'Injector',
        'Chamber + Nozzle',
        'Splitter',
        'Turbine Exhaust',
        'Gas Generator',
        '2nd Fuel Pump',
        'Electric Motor',
        'Inverter',
        'Battery',
        'Battery Cooler',
        'Pressurant Tank',
        'Pressurant',
        'Fuel',
        'Oxidizer',
    )
    for i, (value, name) in enumerate(zip(values,mass_names), start=1):
        # yield f'{i:>3}: {str(value)[:8]:>10} kg - {name}'
        yield f'{i:>3}: {value:>8.1f} kg - {name}'


if __name__ == '__main__':
    from EngineArguments import arguments as args

    design_args = {'thrust': 100e3,
                   'burn_time': 390,
                   'combustion_chamber_pressure': 10e6,
                   'is_frozen': True,
                   'ambient_pressure': None,
                   'exit_pressure_forced': .002e6, }

    cycle_list = (
        (ElectricPumpCycle, args.ep_arguments),
        # (GasGeneratorCycle, args.gg_arguments),
        # (CoolantBleedCycle, args.cb_arguments),
        # (OpenExpanderCycle, args.oe_arguments),
    )

    for Cycle, extra_args in cycle_list:
        complete_args = args.base_arguments | extra_args | design_args | {'verbose': True}
        engine = Cycle(**complete_args)
        make_mass_schematic(engine)
