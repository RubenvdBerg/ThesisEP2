from BaseEngineCycle.OpenCycle import OpenEngineCycle
from dataclasses import dataclass

@dataclass
class ClosedExpandercycle(OpenEngineCycle):
    # TODO: dataclasses inheritance is stupid, see EP-class
    pass_attribute: float = 0
    
    def __post_init__(self):
        
        super().__post_init__()

    @property
    def a(self):
        return 