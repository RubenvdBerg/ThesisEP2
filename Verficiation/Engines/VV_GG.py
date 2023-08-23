if __name__ == '__main__':
    from EngineCycles.GasGeneratorCycle import GasGeneratorCycle, GasGeneratorCycle_DoubleTurbine, \
        GasGeneratorCycle_DoubleTurbineSeries
    from plots.Imaging.performance_image import make_performance_schematic
    from math import sqrt, pi
    import pandas as pd
    import numpy as np
    from time import strftime
    import os
    from EngineArguments.GasGeneratorEngines import hm7b_kwargs, h1_kwargs, rs27_kwargs, f1_kwargs, j2_kwargs, \
        hm60_kwargs
    hm60_kwargs['ambient_pressure'] = None
    engine_dict = {
        'HM7-B': ('single', hm7b_kwargs),
        'H-1': ('single', h1_kwargs),
        'RS-27': ('single', rs27_kwargs),
        'F-1': ('single', f1_kwargs),
        'HM60': ('double', hm60_kwargs),
        'J-2': ('double_single', j2_kwargs),
    }
    engine_class_selector = {'single': GasGeneratorCycle,
                             'double': GasGeneratorCycle_DoubleTurbine,
                             'double_single': GasGeneratorCycle_DoubleTurbineSeries}
    data = {}
    for engine_name, (engine_type, engine_kwargs) in engine_dict.items():
        engine_kwargs.pop('gg_mass_mixture_ratio')
        Engine_Class = engine_class_selector[engine_type]
        engine: GasGeneratorCycle
        engine = Engine_Class(**engine_kwargs)
        main_data = [engine.thrust_chamber.length,
                     sqrt(engine.exit_area / pi) * 2,
                     engine.engine_dry_mass,
                     engine.combustion_chamber.length + engine.nozzle.conv_length,
                     engine.combustion_chamber.radius * 2,
                     engine.thrust_chamber.mass + engine.injector.mass,
                     engine.gas_generator.outlet_mass_flow,
                     engine.gas_generator.mass_mixture_ratio,
                     engine.main_fuel_flow,
                     engine.main_oxidizer_flow,
                     engine.fuel_pump_expected_pressure * 1e-6,
                     engine.oxidizer_pump_expected_pressure * 1e-6, ]

        if engine_type == 'single':
            extra_data = [engine.pumps_power_required * 1e-3,
                          None,]
        else:
            extra_data = [engine.fuel_pumps_power_required * 1e-3,
                          engine.oxidizer_pumps_power_required * 1e-3,]
        total_data = [*main_data, *extra_data, engine.overall_specific_impulse]
        data[engine_name] = total_data

        # make_performance_schematic(engine)
    index = [
        'Thrust Chamber Length [m]',
        'Exit Diameter [m]',
        'Engine Dry Mass [kg]',
        'Chamber Length [m]',
        'Chamber Diameter [m]',
        'Chamber Mass [kg]',
        'GG Mass Flow [kg/s]',
        'GG Mixture Ratio [-]',
        'Fuel Mass Flow [kg/s]',
        'Oxidizer Mass Flow [kg/s]',
        'Fuel Pump Outlet Pressure [MPa]',
        'Oxidizer Pump Outlet Pressure [MPa]',
        '(Fuel) Turbine Power [kW]',
        'Oxidizer Turbine Power [kW]',
        'Specific Impulse [s]',
    ]
    with open(r'C:\Users\rvand\PycharmProjects\ThesisEP2\data\Data\McHugh_Output_Data_GG.csv') as file:
        df1 = pd.read_csv(file, index_col=0, header=0, usecols=list(range(0, len(engine_dict) + 1)))
    df2 = pd.DataFrame(data, index=index)
    df3 = df2.div(df1).apply(lambda x: abs(x - 1) * 100)
    df3 = df3.replace(np.nan, '-')
    time = strftime("%Y%m%d_%H%M%S")
    filename = rf'C:\Users\rvand\PycharmProjects\ThesisEP2\Verficiation\data\VV_GG_results_{time}.xlsx'

    with pd.ExcelWriter(filename, engine='xlsxwriter') as writer:
        df2.to_excel(writer, sheet_name='Simulated')
        df1.to_excel(writer, sheet_name='McHugh_Data')
        df3.to_excel(writer, sheet_name='Difference')
        workbook = writer.book
        worksheet = writer.sheets['Difference']
        worksheet.conditional_format('A1:Z30', {'type': '3_color_scale',
                                                'min_type': 'num',
                                                'mid_type': 'num',
                                                'max_type': 'num',
                                                'min_value': 0,
                                                'mid_value': 15,
                                                'max_value': 30,
                                                'min_color': '#63BE7B',
                                                'mid_color': '#FFEB84',
                                                'max_color': '#F8696B',
                                                })
    os.system(f'start EXCEL.EXE {filename}')
