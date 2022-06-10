from dataclasses import dataclass, field
from BaseEngineCycle.OpenCycle import OpenEngineCycle
from BaseEngineCycle.Turbine import Turbine


@dataclass
class OpenExpandercycle(OpenEngineCycle):
    # TODO: dataclasses inheritance is stupid, see EP- and GG-class
    pass_attribute: float = 0


    # Non-init variables
    turbine_mass_flow: float = field(init=False, default=0)  # [kg/s]

    def __post_init__(self):
        super().__post_init__()


