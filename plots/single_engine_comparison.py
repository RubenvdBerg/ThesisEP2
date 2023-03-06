from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle_DoublePump
from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.Abstract.OpenCycle import OpenEngineCycle
import EngineArguments.arguments as args
from plots.Imaging.mass_image import make_mass_schematic
from plots.Imaging.performance_image import make_performance_schematic
from typing import Sequence
import matplotlib.pyplot as plt

extra_args_switch = {
    ElectricPumpCycle: args.ep_arguments,
    GasGeneratorCycle: args.gg_arguments,
    OpenExpanderCycle_DoublePump: args.oe_arguments,
}


def get_aggregate_mass_dict(engine: EngineCycle) -> dict:
    engine_type = type(engine)
    if issubclass(engine_type, ElectricPumpCycle):
        extra_dict = {'Battery': engine.battery.mass}
    elif issubclass(engine_type, OpenEngineCycle):
        extra_dict = {'Turbine Prop.': engine.turbine_propellant_mass}
    else:
        raise NotImplementedError

    return {'Tanks': engine.tanks_mass} | extra_dict | {
        'Pressurant': engine.pressurant.mass,
        'Feed System': engine.feed_system_mass,
        'Thrust Chamber': engine.total_thrust_chamber_mass,
        'Chamber Prop.': engine.chamber_propellant_mass,
    }


def plot_mass_pie_chart(mass_dict: dict, title_name: str):
    names, values = zip(*mass_dict.items())
    fig, ax = plt.subplots()
    exploded = [.1 for i in range(len(values))]
    ax.pie(values, labels=names, explode=exploded, autopct='%1.2f%%', shadow=True)
    ax.set_title(f'Mass Distribution - {title_name}')
    plt.show()


def plot_mass_pie_chart_comparison(design_args: dict, exclude_cc_prop: bool = True):
    extra_args_switch = {
        ElectricPumpCycle: args.ep_arguments,
        GasGeneratorCycle: args.gg_arguments,
        OpenExpanderCycle_DoublePump: args.oe_arguments,
    }
    name_switch = {
        ElectricPumpCycle: 'EP',
        GasGeneratorCycle: 'GG',
        OpenExpanderCycle_DoublePump: 'OE',
    }

    mass_dicts = {}
    for CycleClass, unique_args in extra_args_switch.items():
        total_args = args.base_arguments | unique_args | design_args
        engine = CycleClass(**total_args)
        make_performance_schematic(engine)
        make_mass_schematic(engine)
        mass_dicts[CycleClass] = get_aggregate_mass_dict(engine)

    if exclude_cc_prop:
        cc_vals = [mass_dict['Chamber Prop.'] for mass_dict in mass_dicts.values()]
        min_cc_val = min(cc_vals)

    fig, axs = plt.subplots(1, 3, figsize=(20, 7.5))
    plt.rcParams['font.size'] = '16'
    for i, (CycleClass, mass_dict) in enumerate(mass_dicts.items()):
        if exclude_cc_prop:
            extra_cc_prop = mass_dict['Chamber Prop.'] - min_cc_val
            if extra_cc_prop:
                mass_dict['Extra CC Prop.'] = extra_cc_prop
            mass_dict.pop('Chamber Prop.')
        names, values = zip(*mass_dict.items())
        exploded = [.1 for i in range(len(values))]
        wedges, texts, autotexts = axs[i].pie(values, explode=exploded, autopct='%1.2f%%', pctdistance=1.3,
                                              counterclock=False)

        axs[i].set_title(name_switch[CycleClass], pad=40, fontsize=20)
    new_names = []
    for name in names:
        if name == 'Turbine Prop.':
            name = 'Battery/' + name
        new_names.append(name)

    fig.legend(wedges, new_names, ncol=3, loc='lower center')
    plt.show()
