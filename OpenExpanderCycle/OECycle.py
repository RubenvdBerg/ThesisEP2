from dataclasses import dataclass, field

import arguments
from BaseEngineCycle.OpenCycle import OpenEngineCycle
from BaseEngineCycle.Turbine import Turbine


@dataclass
class OpenExpandercycle(OpenEngineCycle):
    # Non-init variables
    turbine_mass_flow: float = field(init=False, default=0)  # [kg/s]

    @property  # Override EngineCycle fuelflow -> increase pump requirements -> increase turbine mass flow -> iterate
    def fuel_flow(self):
        return 1 / (self.mass_mixture_ratio + 1) * self.chamber_mass_flow + self.turbine_mass_flow

    @property
    def turbine_inlet_temperature(self):
        return self.cooling_channels.outlet_temperature


if __name__ == '__main__':
    import arguments as args
    engine = OpenExpandercycle(**arguments.le5a_kwargs, **args.open_arguments)
    engine.thrust_chamber.show_contour()
    engine.heat_exchanger.show_heat_flux()
    print(engine.cooling_channels.outlet_temperature)
    print(engine.simple_specific_impulse)
    print(engine.thrust_chamber.length, engine.nozzle.exit_radius*2)