from OpenExpanderCycle.OECycle import OpenExpandercycle
import arguments as args
from Verficiation.VV_Heat_Transfer_Old import convective_heat_transfer_validation


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