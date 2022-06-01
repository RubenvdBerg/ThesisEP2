from BaseEngineCycle.EngineCycle import EngineCycle
from BaseEngineCycle.Cooling import HeatExchanger
from BaseEngineCycle.ThrustChamber import ThrustChamber
from BaseEngineCycle.Nozzle import BellNozzle
from BaseEngineCycle.CombustionChamber import CombustionChamber, Injector
from scipy import constants
import arguments as args


def convective_heat_transfer_validation():
    arguments = args.base_arguments_o.copy()
    del arguments['fuel_name']
    area_ratio = 22
    p_cc = 55e5
    throat_area = .05**2 * constants.pi
    test_engine = EngineCycle(thrust=100e3,
                              combustion_chamber_pressure=p_cc,
                              mass_mixture_ratio=5.6,
                              area_ratio_forced=area_ratio,
                              burn_time=100,
                              is_frozen=False,
                              fuel_name='LH2_NASA',
                              kwak_fix_cycle_type='ep',
                              fuel_pump_specific_power=15e3,
                              oxidizer_pump_specific_power=20e3,
                              **arguments)


    test_chamber = CombustionChamber(
        material_density=1,
        safety_factor=1,
        yield_strength=1,
        throat_area=throat_area,
        combustion_chamber_pressure=p_cc,
        area_ratio_chamber_throat=(80/50)**2,
        characteristic_length=.89)
    test_nozzle = BellNozzle(throat_area=throat_area,
                             area_ratio=22,
                             chamber_radius=test_chamber.radius,
                             conv_chamber_bend_ratio=args.base_arguments['convergent_chamber_bend_ratio'],
                             conv_throat_bend_ratio=args.base_arguments['convergent_throat_bend_ratio'],
                             conv_half_angle=args.base_arguments['convergent_half_angle'],
                             div_throat_half_angle=args.base_arguments['divergent_throat_half_angle'],
                             div_exit_half_angle=args.base_arguments['divergent_exit_half_angle']
                             )
    test_injector = Injector(1,1,1,1,1,False)
    test_thrust_chamber = ThrustChamber(nozzle=test_nozzle,
                                        chamber=test_chamber,
                                        injector=test_injector,
                                        heat_capacity_ratio=test_engine.cc_hot_gas_heat_capacity_ratio)
    test_heat_exchanger = HeatExchanger(thrust_chamber=test_thrust_chamber,
                                        combustion_temperature=test_engine.combustion_temperature,
                                        combustion_chamber_pressure=p_cc,
                                        mass_flow=test_engine.mass_flow,
                                        dynamic_viscosity=test_engine.cc_hot_gas_dynamic_viscosity,
                                        specific_heat_capacity=test_engine.cc_hot_gas_specific_heat_capacity,
                                        hot_gas_emissivity=args.base_arguments['hot_gas_emissivity'],
                                        heat_capacity_ratio=test_engine.cc_hot_gas_heat_capacity_ratio,
                                        maximum_wall_temperature=args.base_arguments['maximum_wall_temperature'],
                                        thrust_chamber_wall_emissivity=.8,
                                        convective_coefficient_mode='Modified Bartz',
                                        prandtl_number=test_engine.cc_hot_gas_prandtl_number)
    test_thrust_chamber.show_contour()
    test_heat_exchanger.show_heat_flux()

    f_string = '.3e'
    print(f'Parameter         Test  \t  Goal')
    print(f'Thr HeatTrsfr:  {test_heat_exchanger.get_convective_heat_flux(0):{f_string}}\t{71e6:.3e}')


if __name__ == '__main__':
    convective_heat_transfer_validation()