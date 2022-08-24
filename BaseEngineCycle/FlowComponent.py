from dataclasses import dataclass, replace
from FlowState import FlowState


@dataclass
class FlowComponent:
    inlet_flow_state: FlowState

    @property
    def inlet_pressure(self):
        return self.inlet_flow_state.pressure

    @property
    def inlet_temperature(self):
        return self.inlet_flow_state.temperature

    @property
    def temperature_change(self):
        return 0

    @property
    def pressure_change(self):
        return 0

    @property
    def outlet_pressure(self):
        return self.inlet_pressure + self.pressure_change

    @property
    def outlet_temperature(self):
        return self.inlet_temperature + self.temperature_change

    @property
    def output_flow_state(self):
        return self.inlet_flow_state.replace(
            pressure=self.outlet_pressure,
            temperature=self.outlet_temperature,
        )