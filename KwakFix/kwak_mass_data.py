import pandas as pd
import os
import EngineArguments.arguments as args
from typing import Optional
from EngineArguments.default_arguments import get_default_kwargs
from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from KwakFix.KwakFixCycles import KwakEngineCycle, KwakFixElectricPumpCycle, KwakFixGasGeneratorCycle


def get_mass_table_data(EngineClass: EngineCycle, design_kwargs: dict, burn_times: list, give_MR_Isp: bool = False,
                        path: Optional[str] = None):
    default_kwargs = get_default_kwargs(EngineClass)
    info = []
    for burn_time in burn_times:
        design_kwargs['burn_time'] = burn_time
        input_kwargs = default_kwargs | design_kwargs
        engine = EngineClass(**input_kwargs)
        if issubclass(EngineClass, ElectricPumpCycle):
            unique = [None,
                      engine.battery.mass,
                      engine.electric_motor.mass,
                      engine.inverter.mass,
                      None, ]
        elif issubclass(EngineClass, GasGeneratorCycle):
            unique = [engine.gg_propellant_mass,
                      None,
                      None,
                      None,
                      engine.gas_generator.mass, ]

        column = [
            engine.chamber_propellant_mass,
            *unique,
            engine.fuel_pump.mass,
            engine.oxidizer_pump.mass,
            engine.oxidizer_tank.mass,
            engine.fuel_tank.mass,
            engine.pressurant_tank.mass,
            engine.pressurant.mass,
            engine.mass_kwak,
        ]
        if give_MR_Isp:
            column += [
                engine.mass_ratio_kwak,
                engine.overall_specific_impulse,
            ]
        info.append(column)
    indices = [
        'CC Prop.',
        'GG Prop.',
        'Battery',
        'Inverter',
        'Electric Motor',
        'Gas Generator',
        'Fuel Pump',
        'Oxidizer Pump',
        'Oxidizer Tank',
        'Fuel Tank',
        'Pressurant Tank',
        'Pressurant',
        'Total',
    ]
    if give_MR_Isp:
        indices += ['Mass Ratio',
                    'Overall Isp']
    data = list(zip(*info))
    columns = list(burn_times)
    info_df = pd.DataFrame(data=data, index=indices, columns=columns)
    if path:
        with pd.ExcelWriter(path) as writer:
            info_df.to_excel(writer)
        os.system(f'start EXCEL.EXE {path}')


def get_aggregate_mass_table_data(EngineClass: EngineCycle, design_kwargs: dict, burn_times: list, give_MR_Isp: bool = False,
                        path: Optional[str] = None):
    default_kwargs = get_default_kwargs(EngineClass)
    info = []
    for burn_time in burn_times:
        design_kwargs['burn_time'] = burn_time
        input_kwargs = default_kwargs | design_kwargs
        engine = EngineClass(**input_kwargs)
        if issubclass(EngineClass, ElectricPumpCycle):
            unique = engine.battery.mass
        elif issubclass(EngineClass, GasGeneratorCycle):
            unique = engine.gg_propellant_mass,

        engine: EngineCycle
        column = [
            engine.chamber_propellant_mass,
            unique,
            engine.feed_system_mass,
            engine.tanks_mass,
            engine.pressurant.mass,
            engine.mass_kwak,
        ]
        if give_MR_Isp:
            column += [
                engine.mass_ratio_kwak,
                engine.overall_specific_impulse,
            ]
        info.append(column)
    indices = [
        'CC Prop.',
        'Power Source',
        'Feed Sys.',
        'Tanks',
        'Pressurant',
        'Total',
    ]
    if give_MR_Isp:
        indices += ['Mass Ratio',
                    'Overall Isp']
    data = list(zip(*info))
    columns = list(burn_times)
    info_df = pd.DataFrame(data=data, index=indices, columns=columns)
    if path:
        with pd.ExcelWriter(path) as writer:
            info_df.to_excel(writer)
        os.system(f'start EXCEL.EXE {path}')

if __name__ == '__main__':
    design_args = {
        'thrust': 100e3,
        'combustion_chamber_pressure': 10e6,
    }
    # get_aggregate_mass_table_data(EngineClass=KwakFixElectricPumpCycle, design_kwargs=design_args, burn_times=[300, 390, 1200],
    #                     give_MR_Isp=True, path='kwak_fix_ep_30-39-120_aggr.xlsx')
    get_mass_table_data(EngineClass=KwakFixElectricPumpCycle, design_kwargs=design_args, burn_times=[300, 390, 1200],
                        give_MR_Isp=True, path='kwak_fix_ep_30-39-120.xlsx')
