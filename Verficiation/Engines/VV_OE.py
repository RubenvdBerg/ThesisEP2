from EngineCycles.OpenExpanderCycle import OpenExpanderCycle_DoubleTurbineSeries
from EngineArguments.ExpanderEngines import le5b_kwargs
from math import pi
import pandas as pd
from time import strftime
import os

engine = OpenExpanderCycle_DoubleTurbineSeries(**le5b_kwargs)
engine_data = {
    'Vacuum Thrust [N]': engine.thrust,
    'O/F-ratio [-]': engine.mass_mixture_ratio,
    'Expansion Ratio [-]': engine.expansion_ratio,
    'Combustion Pressure [Pa]': engine.combustion_chamber_pressure,
    'Burn Time [s]': engine.burn_time,
    'Exit Diameter [m]': (engine.exit_area * 4 / pi) ** .5,
    'Fuel Pump Inlet Temperature [K]': engine.fuel_pump.inlet_temperature,
    'Oxidizer Pump Inlet Temperature [K]': engine.oxidizer_pump.inlet_temperature,
    'Vacuum Specific Impulse [s]': engine.overall_specific_impulse,
    'Length [m]': engine.thrust_chamber.length,
    'Mass [kg]': engine.engine_dry_mass,
    'Fuel Pump Discharge Pressure [Pa]': engine.fuel_pump.outlet_pressure,
    'Oxidizer Pump Discharge Pressure [Pa]': engine.oxidizer_pump.outlet_pressure,
    'Fuel Pump Mass Flow [kg/s]': engine.fuel_pump.mass_flow,
    'Oxidizer Pump Mass Flow [kg/s]': engine.oxidizer_pump.mass_flow,
    'Fuel Shaft Power [W]': engine.fuel_pumps_power_required,
    'Oxidizer Shaft Power [W]': engine.oxidizer_pumps_power_required,
    'Fuel Turbine Inlet Pressure [Pa]': engine.fuel_turbine.inlet_pressure,
    'Oxidizer Turbine Inlet Pressure [Pa]': engine.oxidizer_turbine.inlet_pressure,
    'Fuel Turbine Inlet Temperature [K]': engine.fuel_turbine.inlet_temperature,
    'Oxidizer Turbine Inlet Temperature [K]': engine.oxidizer_turbine.inlet_temperature,
    'Fuel Turbine Outlet Pressure [Pa]': engine.fuel_turbine.outlet_pressure,
    'Oxidizer Turbine Outlet Pressure [Pa]': engine.oxidizer_turbine.outlet_pressure,
    'Fuel Turbine Mass Flow [kg/s]': engine.fuel_turbine.mass_flow,
    'Oxidizer Turbine Mass Flow [kg/s]': engine.oxidizer_turbine.mass_flow_required,
    'Fuel Turbine Pressure Ratio [-]': engine.fuel_turbine.pressure_ratio,
    'Oxidizer Turbine Pressure Ratio [-]': engine.oxidizer_turbine.pressure_ratio,
    'Heat Transfer [W]': engine.heat_flow_rate
}

data_path = r'C:\Users\rvand\PycharmProjects\ThesisEP2\data\Data\Aoki_2001_LE5B_Data.csv'
vv_df = pd.read_csv(data_path, index_col=0, names=['Aoki'])
engine_df = pd.DataFrame.from_dict(engine_data, orient='index', columns=['RoCAT'])
combined_df = vv_df.merge(engine_df, left_index=True, right_index=True, how='outer')
combined_df[r'Diff [\%]'] = abs(combined_df['RoCAT'] / combined_df['Aoki'] - 1) * 100

time = strftime("%Y%m%d_%H%M%S")
filename = rf'C:\Users\rvand\PycharmProjects\ThesisEP2\Verficiation\data\LE5B_results_{time}.xlsx'
with pd.ExcelWriter(filename, engine='xlsxwriter') as writer:
    combined_df.to_excel(writer, header=True, index=True)
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']
    format1 = workbook.add_format({'num_format': '0.00'})
    worksheet.conditional_format('D2:D30', {'type': '3_color_scale',
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
    worksheet.set_column(3, 3, None, format1)
os.system(f'start EXCEL.EXE {filename}')
