from dataclasses import dataclass
from typing import Optional
from EngineComponents.Abstract.FlowState import FlowState

@dataclass
class Propellant:
    initial_flow_state: FlowState
    burn_time: float  # [s]
    margin_factor: float  # [-]
    propellant_density: Optional[float]  # [kg/m3]

    @property
    def name(self):
        return self.initial_flow_state.propellant_name

    @property
    def mass_flow(self):
        return self.initial_flow_state.mass_flow

    @property
    def density(self):
        if self.propellant_density is None:
            return self.initial_flow_state.density
        else:
            return self.propellant_density

    @property
    def mass(self):
        return self.mass_flow * self.burn_time * self.margin_factor

    @property
    def volume(self):
        return self.mass / self.density

    @property
    def volumetric_flow(self):
        return self.mass_flow / self.density



