from dataclasses import dataclass, field
from BaseEngineCycle.OpenCycle import OpenEngineCycle
from BaseEngineCycle.Turbine import Turbine


@dataclass
class OpenExpandercycle(OpenEngineCycle):
    # TODO: dataclasses inheritance is stupid, see EP- and GG-class
    # Non-init variables
    turbine_mass_flow: float = field(init=False, default=0)  # [kg/s]

    def __post_init__(self):
        super().__post_init__()

    @property  # Override EngineCycle fuelflow -> increase pump requirements -> increase turbine mass flow -> iterate
    def fuel_flow(self):
        return 1 / (self.mass_mixture_ratio + 1) * self.chamber_mass_flow + self.turbine_mass_flow

    @property
    def turbine_inlet_temperature(self):
        return self.cooling_channels.outlet_temperature

if __name__ == '__main__':
    import arguments as args
    base_args = args.base_arguments_o
    del base_args['fuel_name']
    del base_args['exit_pressure_forced']
    engine1 = OpenExpandercycle(thrust=121.50E3,
                               combustion_chamber_pressure=40E5,
                               mass_mixture_ratio=5,
                               burn_time=609,
                               fuel_name='LH2_NASA',
                               is_frozen=False,
                               expansion_ratio=130,
                               **base_args, **args.oe_arguments)
    engine2 = OpenExpandercycle(thrust=50E3,
                                combustion_chamber_pressure=10E6,
                                mass_mixture_ratio=2.45,
                                )
    print(engine.cooling_channels.outlet_temperature)
    print(engine.simple_specific_impulse)
    print(engine.thrust_chamber.length, engine.nozzle.exit_radius*2)