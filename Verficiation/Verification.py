from BaseEngineCycle.EngineCycle import EngineCycle
from BaseEngineCycle.HeatExchanger import HeatExchanger
from BaseEngineCycle.ThrustChamber import ThrustChamber
from BaseEngineCycle.Nozzle import BellNozzle
from BaseEngineCycle.CombustionChamber import CombustionChamber, Injector
from OpenExpanderCycle.OECycle import OpenExpandercycle
from scipy import constants
import arguments as args


def convective_heat_transfer_validation():
    arguments = args.base_arguments_o.copy()
    del arguments['fuel_name']
    del arguments['exit_pressure_forced']
    area_ratio = 22
    p_cc = 55e5
    throat_area = .05**2 * constants.pi
    test_engine = EngineCycle(thrust=100e3,
                              combustion_chamber_pressure=p_cc,
                              mass_mixture_ratio=5.6,
                              expansion_ratio=area_ratio,
                              burn_time=100,
                              is_frozen=False,
                              fuel_name='LH2_NASA',
                              fuel_pump_specific_power=15e3,
                              oxidizer_pump_specific_power=20e3,
                              **arguments)
    y_cc, cp_cc = test_engine.cc_hot_gas_heat_capacity_ratio, test_engine.cc_hot_gas_specific_heat_capacity
    mu_cc, pr_cc = test_engine.cc_hot_gas_dynamic_viscosity, test_engine.cc_hot_gas_prandtl_number
    t_c, m_flow = test_engine.combustion_temperature, test_engine.mass_flow
    cea_run = False
    if cea_run:
        y_cc, cp_cc, mu_cc, pr_cc, t_c, m_flow = 1.1441, 8.3713e3, 1.0276e-4, 0.5175, 3395.68, 100e3/4082.3
    test_chamber = CombustionChamber(
        material_density=1,
        safety_factor=1,
        yield_strength=1,
        throat_area=throat_area,
        combustion_chamber_pressure=p_cc,
        # area_ratio_chamber_throat=None,
        area_ratio_chamber_throat=(80 / 50) ** 2,
        propellant_mix='LOX/LH2',
        characteristic_length=None)
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
                                        heat_capacity_ratio=y_cc)
    test_heat_exchanger = HeatExchanger(thrust_chamber=test_thrust_chamber,
                                        combustion_temperature=t_c,
                                        combustion_chamber_pressure=p_cc,
                                        mass_flow=m_flow,
                                        dynamic_viscosity=mu_cc,
                                        specific_heat_capacity=cp_cc,
                                        hot_gas_emissivity=args.base_arguments['hot_gas_emissivity'],
                                        heat_capacity_ratio=y_cc,
                                        maximum_wall_temperature=args.base_arguments['maximum_wall_temperature'],
                                        thrust_chamber_wall_emissivity=.8,
                                        convective_coefficient_mode='Modified Bartz',
                                        prandtl_number=pr_cc)
    test_thrust_chamber.show_contour()
    test_heat_exchanger.show_heat_flux()

    f_string = '.3e'
    print(f'Parameter         Test  \t  Goal')
    print(f'Thr HeatTrsfr:  {test_heat_exchanger.get_convective_heat_flux(0):{f_string}}\t{71e6:.3e}')
    print(test_heat_exchanger.total_heat_transfer)
    print(test_engine.fuel_flow)


def open_expander_verification():
    base_args = args.base_arguments_o.copy()
    del base_args['fuel_name']
    del base_args['exit_pressure_forced']
    engine = OpenExpandercycle(thrust=121.50E3,
                               combustion_chamber_pressure=40E5,
                               mass_mixture_ratio=5,
                               burn_time=609,
                               fuel_name='LH2_NASA',
                               is_frozen=False,
                               expansion_ratio=130,
                               **base_args, **args.oe_arguments)

    print(engine.simple_specific_impulse)
    print(engine.thrust_chamber.length, engine.nozzle.exit_radius*2)


if __name__ == '__main__':
    convective_heat_transfer_validation()
    # open_expander_verification()