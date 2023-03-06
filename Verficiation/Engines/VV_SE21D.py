import math
from dataclasses import dataclass
from EngineArguments import arguments as args
from math import log10, radians
import os
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle_DoublePumpTurbine
from EngineComponents.Abstract.Material import NarloyZ


@dataclass
class SE21D_V3(OpenExpanderCycle_DoublePumpTurbine):
    """See 'Sippel et al. 2003 - Studies on Expander Bleed Cycle Engines for Launchers' for this rocket's configuration.
        """
    @property
    def turbine_mass_flow_initial_guess(self):
        return self.base_mass_flow * .015


@dataclass
class SE21D_Exact_V3(SE21D_V3):
    cooling_section_pressure_drop: float = 3.36e6

    @property
    def fuel_pump_expected_pressure(self):
        return 8.749e6

    @property
    def oxidizer_pump_expected_pressure(self):
        return 8.749e6

    @property
    def secondary_fuel_pump_expected_pressure(self):
        return 12.109e6

    @property
    def throat_area(self):
        return math.pi * .286 ** 2

    def set_heat_transfer(self):
        self.heat_flow_rate = 111.037e6


se_21d_kwargs = args.base_arguments_o | {
    'thrust': 1947.03e3,
    'combustion_chamber_pressure': 6.649e6,
    'fuel_initial_temperature': 21,
    'oxidizer_initial_temperature': 90,
    'expansion_ratio': 12.52,
    'mass_mixture_ratio': 5.5,
    'fuel_name': 'LH2_NASA',
    'burn_time': 100,
    'expansion_ratio_end_cooling': 5,
    'chamber_characteristic_length': 4.0,
    'fuel_pump_efficiency': .7,
    'secondary_fuel_pump_efficiency': .75,
    'oxidizer_pump_efficiency': .76,
    'fuel_initial_pressure': .3e6,
    'oxidizer_initial_pressure': .5e6,
    'turbine_maximum_temperature': 506.452,
    'turbopump_specific_power': 13.5E3,
    'area_ratio_chamber_throat': (.985 / 2) ** 2 / .286 ** 2,
    'ambient_pressure': 101325,
    'divergent_throat_half_angle': radians(15),
    'specific_impulse_quality_factor': .99,
    'oxidizer_turbine_efficiency': .45,
    'fuel_turbine_efficiency': .45,
    'shaft_mechanical_efficiency': .99,
    'oxidizer_exhaust_exit_pressure_forced': .04e6,
    'fuel_exhaust_exit_pressure_forced': .04e6,
    'oxidizer_turbine_outlet_pressure_forced': .3e6,
    'fuel_turbine_outlet_pressure_forced': .3e6,
    'oxidizer_secondary_specific_impulse_quality_factor': .98,
    'fuel_secondary_specific_impulse_quality_factor': .98,
    'exhaust_material': NarloyZ,
}

se_21d_vac_kwargs = se_21d_kwargs | {'thrust': 2196.682e3, 'ambient_pressure': 0}


def get_SE21D_data(is_pressure_exact: bool, is_vacuum: bool = False, show_schematic: bool = False):
    engineclass = SE21D_Exact_V3 if is_pressure_exact else SE21D_V3
    base_args = se_21d_vac_kwargs if is_vacuum else se_21d_kwargs
    engine = engineclass(**base_args)

    if show_schematic:
        from plots.Imaging.performance_image import make_performance_schematic
        engine.engine_dry_mass
        make_performance_schematic(engine)

    name = '2nd Estimate' if is_pressure_exact else 'Estimate'
    return (('Name', 'Unit', name, 'Expected'),
            ('Chamber Sp. Impulse', '[s]', engine.chamber_specific_impulse, 365.085),
            ('Exhaust Sp. Impulse Fu.', '[s]', engine.oxidizer_secondary_exhaust.specific_impulse, 158.933),
            ('Exhaust Sp. Impulse Ox.', '[s]', engine.fuel_secondary_exhaust.specific_impulse, 153.064),
            ('Overall Sp. Impulse', '[s]', engine.overall_specific_impulse, 360.985),
            (r'$\Delta P$ Oxid. Pump', '[MPa]', engine.oxidizer_pump.pressure_change * 1e-6, 8.749 - .5),
            (r'$\Delta P$ Fuel Pump', '[MPa]', engine.fuel_pump.pressure_change * 1e-6, 8.749 - .3),
            (r'$\Delta P$ Fuel Pump 2', '[MPa]', engine.secondary_fuel_pump.pressure_change * 1e-6, 12.109 - 8.749),
            ('Fuel Pump Power Req.', '[MW]', engine.fuel_pump.power_required * 1e-6, 15.443),
            ('Fuel Pump 2 Power Req.', '[MW]', engine.secondary_fuel_pump.power_required * 1e-6, 1.01),
            ('Oxid. Pump Power Req.', '[MW]', engine.oxidizer_pump.power_required * 1e-6, 4.315),
            ('Total Pump Power Req.', '[MW]', engine.pumps_power_required * 1e-6, 20.768),
            # ('TC Length', '[m]', engine.thrust_chamber.length, 0),
            ('Heat Transfer', '[MW]', engine.heat_flow_rate * 1e-6, 111.037),
            # ('Turb. In. Temp.', '[K]', engine.turbine_inlet_temperature, 506.452),
            ('Turb. Mass Flow', '[kg/s]', engine.turbine_mass_flow, 10.174),
            ('Cool. Mass Flow', '[kg/s]', engine.cooling_channel_section.inlet_mass_flow, 16.394),
            ('Main Fuel Flow', '[kg/s]', engine.main_fuel_flow, 93.677),
            ('Main Oxid. Flow', '[kg/s]', engine.main_oxidizer_flow, 456.323),
            ('Chamber Diameter', '[m]', engine.combustion_chamber.radius * 2, 0.985),
            ('Chamber Volume', '[m$^3$]', engine.combustion_chamber.volume_incl_nozzle_convergent, 1.029),
            ('Subs. Length', '[m]', engine.combustion_chamber.length + engine.nozzle.conv_length, 1.582),
            ('Throat Radius', '[m]', engine.nozzle.throat_radius, 0.286),
            ('Nozzle Length', '[m]', engine.nozzle.div_length, 2.548),
            )


def print_SE21_data(data: tuple):
    for name, unit, value, expected in data:
        print(f'{name:<18} \t {unit:<7} \t {value:<8.3f} \t {expected:.3f}')


def format_to_n_digits(number: float, n: int = 5):
    try:
        n_before_comma = int(log10(abs(number)) // 1 + 1) if number != 0 else 0
        decimals = n - n_before_comma
        return f'{number:.{decimals}f}'.lstrip('0')
    except TypeError:
        return number


def write_SE21D_data_to_excel(**kwargs):
    data_estimated = get_SE21D_data(is_pressure_exact=False, **kwargs)
    data_exact = get_SE21D_data(is_pressure_exact=True, **kwargs)

    def diff(new_list: list, old_list: list):
        return [r'Diff. [\%]'] + [None if new is None else abs(new - old) / old * 100 for new, old in
                                zip(new_list[1:], old_list[1:])]

    var_col = [row[0] for row in data_estimated]  # Variable
    unit_col = [row[1] for row in data_estimated]  # Unit
    exp_col = [row[3] for row in data_estimated]  # Expected
    empty_col = [None] * 20
    exact_col = [row[2] for row in data_exact]  # Exact
    d_exact_col = diff(exact_col, exp_col)  # Exact Difference
    estm_col = [row[2] for row in data_estimated]  # Estimated
    d_estm_col = diff(estm_col, exp_col)  # Estimated Difference

    import pandas as pd
    from itertools import zip_longest
    from time import strftime

    exact_col = [format_to_n_digits(x) for x in exact_col]
    estm_col = [format_to_n_digits(x) for x in estm_col]

    columns = zip_longest(
        var_col, unit_col, exp_col, estm_col, d_estm_col, exact_col, d_exact_col, fillvalue=None
    )

    df = pd.DataFrame(data=columns)
    time = strftime("%Y%m%d_%H%M%S")
    filename = rf'C:\Users\rvand\PycharmProjects\ThesisEP2\Verficiation\data\SE21D_results_{time}.xlsx'
    with pd.ExcelWriter(filename, engine='xlsxwriter') as writer:
        df.to_excel(writer, header=False, index=False)
        workbook = writer.book
        worksheet = writer.sheets['Sheet1']
        format1 = workbook.add_format({'num_format': '0.00'})
        format2 = workbook.add_format({'num_format': '0.000'})
        for range in ['G2:G5', 'G9:G12', 'G14:G20', 'G22', 'E2:E30']:
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

        # worksheet.set_column(0, 8, None, format2)
        worksheet.set_column(4, 4, None, format1)
        worksheet.set_column(6, 6, None, format1)
    os.system(f'start EXCEL.EXE {filename}')


if __name__ == '__main__':
    write_SE21D_data_to_excel(show_schematic=True, is_vacuum=False)


