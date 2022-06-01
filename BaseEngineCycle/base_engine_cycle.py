from math import pi
from BaseEngineCycle.EngineCycle import EngineCycle

if __name__ == '__main__':
    from arguments import base_arguments_o

    throat_area = pi * .05 ** 2
    chamber_pressure = 55e5
    area_ratio = 22
    del base_arguments_o['exit_pressure']
    b = EngineCycle(thrust=100e3,
                    burn_time=160,
                    combustion_chamber_pressure=chamber_pressure,
                    mass_mixture_ratio=5.6,
                    is_frozen=False,
                    fuel_pump_specific_power=15E3,
                    oxidizer_pump_specific_power=20E3,
                    kwak_fix_cycle_type='ep',
                    **base_arguments_o,
                    exit_pressure=.5e5,
                    )

    # b.thrust_chamber.show_contour()
    # b.thrust_chamber.show_mach()
    # b.heat_exchanger.show_adiabatic_wall_temp()
    # b.heat_exchanger.show_heat_flux()
    # b.heat_exchanger.show_heat_flux_coefficient()
    # print(b.thrust_chamber.surface)
    print(f'Conv. transfer: {b.heat_exchanger.total_convective_heat_transfer:.3E}')
    print(f'Radi. transfer: {b.heat_exchanger.total_radiative_heat_transfer:.3E}')
    print(f'Total transfer: {b.heat_exchanger.total_heat_transfer}')
    # 0.7284979653537342
    # Heat
    # transfer: 2.822E+07