from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.SimpleElectricPumpCycle import SimpleElectricPumpCycle
from EngineArguments.ElectricPumpEngines import rutherford_kwargs, lee_kwargs
from plots.Imaging.performance_image import make_performance_schematic
from plots.Imaging.mass_image import make_mass_schematic
import pandas as pd
import os
from time import strftime

engine = SimpleElectricPumpCycle(**lee_kwargs)


isp = 500/ (97.7883/1.1/600*9.80665)
print(isp, engine.overall_specific_impulse, f'{isp/engine.overall_specific_impulse-1:.3%}')
lee_mass_values = {
    'Tanks': .5977,
    'Feed line': 2.4735,
    'Power': 2.5853,
    'Thrust Chamber': 5.1382,
    'Propellant': 97.7883,

}

rocat_mass_values = {
    'Tanks': engine.tanks_mass,
    'Feed line': 0,
    'Power': engine.battery.mass + engine.feed_system_mass,
    'Thrust Chamber': engine.thrust_chamber.mass + engine.injector.mass,
    'Propellant': engine.props_mass,
}

for dict in [rocat_mass_values, lee_mass_values]:
    dict['Total'] = sum(dict.values()) - dict['Feed line']
    dict['Dry'] = dict['Total'] - dict['Propellant']


def diff(val1, val2):
    return abs(val2 / val1 - 1) * 100


data = {key: [value, rocat_mass_values[key], diff(value, rocat_mass_values[key])] for key, value in
        lee_mass_values.items()}
df = pd.DataFrame.from_dict(data=data, orient='index')
# for key, lee_value in lee_mass_values.items():
#     rocat_value = rocat_mass_values[key]
#     print(f'{key:<15} - Lee:{lee_value:>10.4f}, RoCAT:{rocat_value:>10.4f}, Diff:{rocat_value / lee_value - 1:>8.2%}')
# print(engine.initial_mass, engine.dry_mass)
time = strftime("%Y%m%d_%H%M%S")
filename = rf'C:\Users\rvand\PycharmProjects\ThesisEP2\Verficiation\data\VV_EP_Lee_results_{time}.xlsx'
with pd.ExcelWriter(filename, engine='xlsxwriter') as writer:
    df.to_excel(writer, header=['Lee', 'RoCAT', r'Diff. [\%]'], index=True)
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']
    format1 = workbook.add_format({'num_format': '0.00'})
    format2 = workbook.add_format({'num_format': '0.000'})
    for range in ['D2:D8',]:
        worksheet.conditional_format(range, {'type': '3_color_scale',
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
    worksheet.set_column(1, 2, None, format2)
    worksheet.set_column(3, 3, None, format1)
os.system(f'start EXCEL.EXE {filename}')
